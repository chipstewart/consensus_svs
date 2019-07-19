function [ X ] = MantaVCF2dRangerForBP(VCF,pair_id, normpaneldb,use_PoN_JH,minscoreforbp,build,MINSPAN,MAXPON,MAXNORM,MINTUMSR,MINTUMRP,MINTUM)
%function [ X ] = MantaVCF2dRangerForBP(VCF,pair_id,normpaneldb,use_PoN_JH,minscoreforbp,build,MINSPAN,MAXPON,MAXNORM,MINTUMSR,MINTUMRP,MINTUM) 
%   Convert Manta SV  VCF to a dRanger forBP tsv 
%   as close as  possible.
%   CS 05 Dec 2017
%
VCF
pair_id
normpaneldb
use_PoN_JH

if ischar(VCF)
    cmd=['cp -f ' VCF '  .']; unix(cmd)
    [p1 f1 e1 ]=fileparts(VCF)
    if ismember(e1,'.gz')
        cmd=['gunzip -f  ./' f1  e1]; unix(cmd)
    else
        f1=[f1 e1];
    end
    V=load_vcf(['./' f1])    
elseif isstruct(VCF)
    V=VCF;
else
    return
end
individual=pair_id;

if ~exist('MINSPAN','var'), MINSPAN = 200; end
if ~exist('MAXPON','var'), MAXPON = 0; end
if ~exist('MAXNORM','var'), MAXNORM = 0; end
if ~exist('MINTUM','var'), MINTUM = 4; end
if ~exist('MINTUMRP','var'), MINTUMRP = 1; end
if ~exist('MINTUMSR','var'), MINTUMSR = 1; end
if ~exist('build','var'), build = 'hg19'; end
if ~exist('normpaneldb','var'), normpaneldb = ''; end
if ~exist('minscoreforbp','var'), minscoreforbp = 0.01; end
if ~exist('use_PoN_JH','var'), use_PoN_JH = 'TRUE'; end
if ~exist('min_ignore_matchingnormal','var'), min_ignore_matchingnormal = Inf; end

if ~isnumeric(MINSPAN), MINSPAN = str2double(MINSPAN); end
if ~isnumeric(MAXPON), MAXPON = str2double(MAXPON); end
if ~isnumeric(MAXNORM), MAXNORM = str2double(MAXNORM); end
if ~isnumeric(MINTUM), MINTUM = str2double(MINTUM); end
if ~isnumeric(MINTUMRP), MINTUMRP = str2double(MINTUMRP); end
if ~isnumeric(MINTUMSR), MINTUMSR = str2double(MINTUMSR); end

if ~isnumeric(minscoreforbp), minscoreforbp = str2double(minscoreforbp); end
if ischar(use_PoN_JH), use_PoN_JH=ismember(upper(cellstr(use_PoN_JH)),'TRUE')>0; end

P=[];
P = impose_default_value(P,'dRanger_input_iszfile','all.isz');
P = impose_default_value(P,'panel_of_normals_database',normpaneldb);
P = impose_default_value(P,'results_name','manta_results');
P = impose_default_value(P,'min_score_to_perform_BreakPointer',minscoreforbp);
P = impose_default_value(P,'build',build);
P = impose_default_value(P,'use_new_method',use_PoN_JH);
P = impose_default_value(P,'build_refseq',P.build);
P = impose_default_value(P,'impute_promoters',true);
P = impose_default_value(P,'imputed_promoter_size',3000);
P = impose_default_value(P,'min_ignore_matchingnormal',min_ignore_matchingnormal);

if P.use_new_method
    P.panel_of_normals_directory=P.panel_of_normals_database;
end


if length(V.CHROM)<1
    L='individual	num	chr1	str1	pos1	chr2	str2	pos2	class	span	tumreads	normreads	normpanelbins	min1max1	range1	stdev1	min2	max2	range2	stdev2	gene1	site1	gene2	site2	fusion	fmapqzT1	fmapqzN1	fmapqzT2	fmapqzN2	nuwpT1	nuwpN1	nuwpT2	nuwpN2	zstdev1	zstdev2	quality	score	somatic	somatic_score	BPtry	t_alt_RP	t_alt_SR	t_depn_alt_RP	n_alt_SR	n_dep'
    fname = [individual '.' P.results_name '.forBP.txt'];
    fid=fopen(fname,'wt');
    fprintf(fid,'%s\n',L);
    fclose(fid);
    X=[];
    return
end
    
% hard-codeed parameters
P.minpairs = 0;
P.minsomratio = 10;

disp(P)

ff=fieldnames(V);
ff=ff(end+[-1:0]);
V.DATA_N=V.(ff{1});
V.DATA_T=V.(ff{2});
V=rmfield(V,ff);
ff=regexprep(ff,'^v__','');
V.DATA_N=strcat(V.DATA_N,':',ff{1});
V.DATA_T=strcat(V.DATA_T,':',ff{2});
V.pair_id = repmat({pair_id},size(V.CHROM));
if isnumeric(V.CHROM)
    V.CHROM=regexprep(cellstr(num2str(V.CHROM)),' ','');
end    

% remove anything outside of 1-22,x,y 
C=[regexprep(cellstr(num2str((1:22)')),' ','');{'X'};{'Y'}]
V=trimStruct(V,ismember(V.CHROM,C))



VE=trimStruct(V,strfindk(V.ALT,'<'));
VX=trimStruct(V,strfindk(V.ALT,'<','v'));
VX=trimStruct(VX,strfindk(VX.INFO,'SVTYPE=BND','v'));

V=trimStruct(V,strfindk(V.INFO,'SVTYPE=BND'));
%V=trimStruct(V,strfindk(V.ALT,'<','v'));
%V.num=str2double(cellfun(@(x) x(2),regexp(V.ID,':','split')));
V.side=1+str2double(cellfun(@(x) x(end),regexp(V.ID,':','split')));
[~,k]=sort(V.ID);
V=trimStruct(V,k);
V.ID0=regexprep(V.ID,':.$','');
t=tab(V.ID0)
q=regexp(V.INFO,'MATEID=.+\;','match')
q=cellfun(@(x) x{1}, regexp([q{:}]',';','split'),'UniformOutput', false)
q=regexprep(q,'MATEID=','')
V.MID=q;

[i1 m1]=ismember(V.ID,V.MID);
V=trimStruct(V,i1)

kpp=cellfun(@length,regexp(V.ALT,'.+].+]','match'))>0;
kpm=cellfun(@length,regexp(V.ALT,'.+[.+[','match'))>0;
kmp=cellfun(@length,regexp(V.ALT,'].+].+','match'))>0;
kmm=cellfun(@length,regexp(V.ALT,'[.+[.+','match'))>0;
V.str1=zeros(size(V.POS));
V.str2=zeros(size(V.POS));
% + str = 0,  - str = 1; 
V.str1(find(kmp))=1;
V.str1(find(kmm))=1;
V.str2(find(kpm))=1;
V.str2(find(kmm))=1;

%q=regexp(V.INFO,'STRANDS=..:','match'); q= [q{:}]'; q=regexprep(q,'STRANDS=',''); q=regexprep(q,':','');
%V.STRANDS=q;

if isnumeric(V.CHROM)
    V.CHROM=num2str(V.CHROM);
end
V1=trimStruct(V,V.side==1);
V2=trimStruct(V,V.side==2);

X=[];

% convert V to forBP format
X.individual=V1.pair_id;
X.num=(1:V1.N)';
X.chr1=chrom2num(V1.CHROM);
X.str1=V1.str1;
X.pos1=V1.POS;
X.chr2=chrom2num(V2.CHROM);
X.str2=V2.str1;
X.pos2=V2.POS;
X.class=repmat({'deletion'},size(X.chr1));
X.span=NaN(size(X.chr1));
k=find(X.chr1==X.chr2);
X.span(k)=X.pos2(k)-X.pos1(k);
k=find( (X.chr1==X.chr2) & X.span<1e6 & X.str2==X.str1  );
X.class(k)={'inversion'};
k=find( (X.chr1==X.chr2) & X.span<1e6 & X.str2==0 &  X.str1==1  );
X.class(k)={'tandem_dup'};
k=find( X.chr1~=X.chr2) ;
X.class(k)={'inter_chr'};
k=find( (X.chr1==X.chr2) & X.span>1e6  );
X.class(k)={'long_range'};


qFMT=regexp(V1.FORMAT,':','split');%qFMT=[qFMT{:}];
qt=regexp(V1.DATA_T,':','split');%qt=[qt{:}];
qn=regexp(V1.DATA_N,':','split');%qn=[qn{:}];
talt=zeros(V1.N,1);
tsr=talt;nalt=talt;nsr=talt;tdep=talt; ndep=talt;
for i=1:V1.N
    k=find(ismember(qFMT{i},'PR'));
    if ~isempty(k)
        pr=qt{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        talt(i)=str2double(pr1(2));
        tdep(i)=str2double(pr1(1));
        pr=qn{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        nalt(i)=str2double(pr1(2));
        ndep(i)=str2double(pr1(1));
    end
    k=find(ismember(qFMT{i},'SR'));
    if ~isempty(k)
        pr=qt{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        tsr(i)=str2double(pr1(2));
        tdep(i)=tdep(i)+str2double(pr1(1));
        pr=qn{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        nsr(i)=str2double(pr1(2));
        ndep(i)=ndep(i)+str2double(pr1(1));
    end
end


X.tumreads=talt;
X.normreads=nalt;
X.normpanelbins=0*X.tumreads;

X.min1=X.pos1-300;
X.max1=X.pos1+300;
X.range1=X.max1-X.min1;
X.stdev1=X.range1/5;

X.min2=X.pos2-300;
X.max2=X.pos2+300;
X.range2=X.max2-X.min2;
X.stdev2=X.range2/5;

X.gene1=repmat({''},size(X.chr1));
X.site1=repmat({''},size(X.chr1));
X.gene2=repmat({''},size(X.chr1));
X.site2=repmat({''},size(X.chr1));
X.fusion=repmat({''},size(X.chr1));


X.fmapqzT1=NaN(size(X.chr1));
X.fmapqzN1=NaN(size(X.chr1));
X.fmapqzT2=NaN(size(X.chr1));
X.fmapqzN2=NaN(size(X.chr1));
X.nuwpT1=NaN(size(X.chr1));
X.nuwpN1=NaN(size(X.chr1));
X.nuwpT2=NaN(size(X.chr1));
X.nuwpN2=NaN(size(X.chr1));
X.zstdev1=NaN(size(X.chr1));
X.zstdev2=NaN(size(X.chr1));
X.quality=NaN(size(X.chr1));
X.score=NaN(size(X.chr1));
X.somatic=ones(size(X.chr1));
X.somatic_score=NaN(size(X.chr1));
X.BPtry=ones(size(X.chr1));

X.quality=ones(size(X.quality));
X.score=X.tumreads.*X.quality;
X.somatic_score=X.score;

X.t_alt_RP=talt;
X.t_alt_SR=tsr;
X.t_dep=tdep;
X.n_alt_RP=nalt;
X.n_alt_SR=nsr;
X.n_dep=ndep;

% event VCF part 
%q=regexp(VE.INFO,'STRANDS=..:','match'); q= [q{:}]'; q=regexprep(q,'STRANDS=',''); q=regexprep(q,':','');
%VE.STRANDS=q;
if VE.N>0
    VE.STRANDS=repmat({'+-'},size(VE.ALT));
    k=find(ismember(VE.ALT,'<DUP:TANDEM>'));
    VE.STRANDS(k)={'-+'}
    %inv=find(ismember(VE.ALT,'<INV>')
    kinv3=find(cellfun(@length,regexp(VE.INFO,'INV3','match'))>0);
    kinv5=find(cellfun(@length,regexp(VE.INFO,'INV5','match'))>0);
    VE.STRANDS(kinv3)={'++'}
    VE.STRANDS(kinv5)={'--'}
    q=regexp(VE.INFO,'END=\w+;','match'); q= [q{:}]'; q=regexprep(q,'END=',''); q=regexprep(q,';','');
    VE.END=str2double(q);
else
    VE.END=VE.POS
    VE.STRANDS=VE.ID
end

X2=trimStruct(X,[]);
X2.num=V1.N+(1:VE.N)';
X2.individual=VE.pair_id;
X2.chr1=chrom2num(VE.CHROM);
X2.str1=zeros(size(X2.chr1));
X2.pos1=VE.POS;
X2.chr2=chrom2num(VE.CHROM);
X2.str2=zeros(size(X2.chr2));
X2.pos2=VE.END;
if VE.N>0
    X2.str1=ismember(substring(VE.STRANDS,1,1),'-');
    X2.str2=ismember(substring(VE.STRANDS,2,1),'-');
end
    
X2.class=repmat({'deletion'},size(X2.chr1));
X2.span=NaN(size(X2.chr1));
k=find(X2.chr1==X2.chr2);
X2.span(k)=X2.pos2(k)-X2.pos1(k);
k=find( (X2.chr1==X2.chr2) & X2.span<1e6 & X2.str2==X2.str1  );
X2.class(k)={'inversion'};
k=find( (X2.chr1==X2.chr2) & X2.span<1e6 & X2.str2==0 &  X2.str1==1  );
X2.class(k)={'tandem_dup'};
k=find( X2.chr1~=X2.chr2) ;
X2.class(k)={'inter_chr'};
k=find( (X2.chr1==X2.chr2) & X2.span>1e6  );
X2.class(k)={'long_range'};


qFMT=regexp(VE.FORMAT,':','split');%qFMT=[qFMT{:}];
qt=regexp(VE.DATA_T,':','split');%qt=[qt{:}];
qn=regexp(VE.DATA_N,':','split');%qn=[qn{:}];
talt=zeros(VE.N,1);
tsr=talt;nalt=talt;nsr=talt;tdep=talt; ndep=talt;
for i=1:VE.N
    k=find(ismember(qFMT{i},'PR'));
    if ~isempty(k)
        pr=qt{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        talt(i)=str2double(pr1(2));
        tdep(i)=str2double(pr1(1));
        pr=qn{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        nalt(i)=str2double(pr1(2));
        ndep(i)=str2double(pr1(1));
    end
    k=find(ismember(qFMT{i},'SR'));
    if ~isempty(k)
        pr=qt{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        tsr(i)=str2double(pr1(2));
        tdep(i)=tdep(i)+str2double(pr1(1));
        pr=qn{i}; pr=pr(k);
        pr1=regexp(pr,',','split');pr1=pr1{1};
        nsr(i)=str2double(pr1(2));
        ndep(i)=ndep(i)+str2double(pr1(1));
    end
end



X2.tumreads=talt;
X2.normreads=nalt;
X2.normpanelbins=0*X2.tumreads;

X2.min1=X2.pos1-300;
X2.max1=X2.pos1+300;
X2.range1=X2.max1-X2.min1;
X2.stdev1=X2.range1/5;

X2.min2=X2.pos2-300;
X2.max2=X2.pos2+300;
X2.range2=X2.max2-X2.min2;
X2.stdev2=X2.range2/5;

X2.gene1=repmat({''},size(X2.chr1));
X2.site1=repmat({''},size(X2.chr1));
X2.gene2=repmat({''},size(X2.chr1));
X2.site2=repmat({''},size(X2.chr1));
X2.fusion=repmat({''},size(X2.chr1));


X2.fmapqzT1=NaN(size(X2.chr1));
X2.fmapqzN1=NaN(size(X2.chr1));
X2.fmapqzT2=NaN(size(X2.chr1));
X2.fmapqzN2=NaN(size(X2.chr1));
X2.nuwpT1=NaN(size(X2.chr1));
X2.nuwpN1=NaN(size(X2.chr1));
X2.nuwpT2=NaN(size(X2.chr1));
X2.nuwpN2=NaN(size(X2.chr1));
X2.zstdev1=NaN(size(X2.chr1));
X2.zstdev2=NaN(size(X2.chr1));
X2.quality=NaN(size(X2.chr1));
X2.score=NaN(size(X2.chr1));
X2.somatic=ones(size(X2.chr1));
X2.somatic_score=NaN(size(X2.chr1));
X2.BPtry=ones(size(X2.chr1));

X2.quality=ones(size(X2.quality));
if VE.N>0
    X2.score=X2.tumreads.*X2.quality;
end
X2.somatic_score=X2.score;

X2.t_alt_RP=talt;
X2.t_alt_SR=tsr;
X2.t_dep=tdep;
X2.n_alt_RP=nalt;
X2.n_alt_SR=nsr;
X2.n_dep=ndep;
X2.STRANDS=VE.STRANDS;
 

X=mergeStruct(X,X2);
X.num=(1:X.N)';

fname1 = [individual '.' P.results_name '.forBP.raw.tsv'];
fprintf('Saving %s\n',fname1);
%save_struct(Xbp,fname);
printStruct(X,-1, fname1);

fprintf('filter low span events >= %d\n',MINSPAN) 
KEEP=sum(~(X.span<MINSPAN))
TOSS=sum(X.span<MINSPAN)

X=trimStruct(X,~(X.span<MINSPAN));

fprintf('normreads <= %d\n',MAXNORM) 
KEEP=sum(X.normreads<=MAXNORM)
TOSS=sum(X.normreads>MAXNORM)
X=trimStruct(X,X.normreads<=MAXNORM);

fprintf('tumor RP + SR reads >= %d\n',MINTUM) 
KEEP=sum((X.t_alt_SR+X.t_alt_RP)>=MINTUM)
TOSS=sum((X.t_alt_SR+X.t_alt_RP)<MINTUM)
X=trimStruct(X,(X.t_alt_SR+X.t_alt_RP)>=MINTUM);


fprintf('tumor RP reads >= %d\n',MINTUMRP) 
KEEP=sum(X.t_alt_RP>=MINTUMRP)
TOSS=sum(X.t_alt_RP<MINTUMRP)
X=trimStruct(X,X.t_alt_RP>=MINTUMRP);

fprintf('tumor SR reads >= %d\n',MINTUMSR) 
KEEP=sum(X.t_alt_SR>=MINTUMSR)
TOSS=sum(X.t_alt_SR<MINTUMSR)
X=trimStruct(X,X.t_alt_SR>=MINTUMSR);


X=rmfield(X,'N')
% ADD CLASS
if length(X.pos1)>0
    X.class = classify_rearrangements(X);
else
    X.class=X.chr1;
end

% ORDER FIELDS
X.individual = repmat({individual},slength(X),1);
flds = {'individual','num','chr1','str1','pos1','chr2','str2','pos2','class','span','tumreads','normreads'};
X = orderfields_first(X,flds);

% SCREEN AGAINST PANEL OF NORMALS (if specified)
if ~isempty(P.panel_of_normals_database)
    if ~isempty(X.num)
        X = dRanger_screen_against_panel_of_normals(X,P);
    else
        X.normpanelbins=X.num;
    end
    if (P.use_new_method)
        flds = [flds 'normpanelbins'];
    else
        flds = [flds 'normpanelreads','normpanelsamps','normpaneldetails'];
    end
    X = orderfields_first(X,flds);
end

% ANNOTATE SITES BY GENE AND FUSION STATUS
X = dRanger_annotate_sites(X,P);

% compute quality
X.quality = dRanger_calculate_quality(X);

% compute score
X.score = X.tumreads .* X.quality;
X.score(round(X.score)<P.minpairs) = 0;     % impose threshold = minpairs

% compute somatic score
if isfield(X,'normpanelbins')
    normreads = max(X.normreads,X.normpanelbins);
    ignore_matchingnormal = (X.tumreads ./ X.normreads) >= P.min_ignore_matchingnormal;
    if any(ignore_matchingnormal)
        normreads(ignore_matchingnormal) = X.normpanelbins(ignore_matchingnormal);
    end
elseif isfield(X,'normpanelreads')
    normreads = max(X.normreads,X.normpanelreads);
    ignore_matchingnormal = (X.tumreads ./ X.normreads) >= P.min_ignore_matchingnormal;
    if any(ignore_matchingnormal)
        normreads(ignore_matchingnormal) = X.normpanelreads(ignore_matchingnormal);
    end
else
    normreads = X.normreads;
end
X.somatic = (X.tumreads ./ normreads) >= P.minsomratio;
X.somatic_score = X.score .* X.somatic;

% sort rearrangements by somatic score
a = nan(slength(X),1);
[X ord] = sort_struct(X,'somatic_score',-1);

% coarse-filter results for BreakPointer
X.BPtry = (X.score>=P.min_score_to_perform_BreakPointer);
Xbp = reorder_struct(X,X.BPtry);

fprintf('normpanelbins<=%d\n',MAXPON) 
KEEP=sum(Xbp.normpanelbins<=MAXPON)
TOSS=sum(Xbp.normpanelbins>MAXPON)
Xbp=trimStruct(Xbp,Xbp.normpanelbins<MAXPON);

% SAVE list of rearrangements for BreakPointer to assemble
fname = [individual '.' P.results_name '.forBP.txt'];
fprintf('Saving %s\n',fname);
%save_struct(Xbp,fname);
printStruct(Xbp,-1, fname);
if isdeployed %(FIREHOSE)
    exit
end
end

function test1
clear
cd ~/Downloads
addpath /Users/stewart/CancerGenomeAnalysis/trunk/matlab
addpath /Users/stewart/CancerGenomeAnalysis/trunk/matlab/seq
addpath /Users/stewart/CancerGenomeAnalysis/trunk/matlab/mike
normpaneldb='/cga/fh/pcawg_pipeline4/refdata/protected/dRanger_PoN_JH'
use_PoN_JH=true
minscoreforbp=0.01
build='/xchip/cga/reference/annotation/db/ucsc/hg19/R.mat'
pair_id='THCA-EM-A3FL-TP-NB'
VCF='/Users/stewart/Downloads/THCA-EM-A3FL-TP-NB.manta.somaticSV.vcf.gz'
MINSPAN=200,MAXPON=1,MAXNORM=2,MINTUMSR=0,MINTUMRP=2,MINTUM=4

X=MantaVCF2dRangerForBP(VCF,pair_id,normpaneldb,use_PoN_JH,minscoreforbp,build,MINSPAN,MAXPON,MAXNORM,MINTUMSR,MINTUMRP,MINTUM)
% V=load_vcf('/Users/stewart/Downloads/lumpy.somatic.sv.vcf');
end
