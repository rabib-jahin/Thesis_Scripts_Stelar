function [rooted_newick,mad_stats]=mad(unrooted_newick,flags)
%
%Minimal Ancestor Deviation Rooting.
%
%Usage:
%
%   [rooted_newick,mad_stats]=mad(unrooted_newick,flags)
%
%   unrooted_newick: Unrooted tree string in NEWICK format.
%   flags          : (optional) string with any of:
%              's' : Append MAD statistics to NEWICK string.
%              't' : Do not round to zero branches < 10^-6.
%
%   rooted_newick: Rooted tree string in NEWICK format.
%   mad_stats    : Rooting statistics and messages.
%
%Version: 2.0, 14-Jan-2018 (stand-alone for distribution)
%


if nargin==0
    z=dbstack;
    help(z(1).name)
    return
end
if nargin<2
    flags='z';
end
rooted_newick='Error';
mad_stats.messages={};
mad_info.flags=flags;
%input
if ~ischar(unrooted_newick) 
    unrooted_newick=[unrooted_newick{:}];
end
if ~ischar(unrooted_newick)
    error('Expecting a NEWICK string.')
end
unrooted_newick=regexprep(unrooted_newick,'\s+','');
n5=arrayfun(@(c)sum(unrooted_newick==c),'(),:;');
if n5(5)~=1 || unrooted_newick(end)~=';'
    error('Corrupt NEWICK string - expeting one '';'' terminated string',s);
end    
if n5(1)~=n5(2)
    error('Corrupt NEWICK string - unbalanced ()',s);
end    
if n5(1)>n5(3)
    error('Corrupt NEWICK string - too many ''()'' / too few '',''.',s);
end    
if n5(4)>n5(3)*2
    error('Corrupt NEWICK string - too many '':'' / too few '',''.',s);
end
tiny=any(upper(flags)=='T');
astt=any(upper(flags)=='S');
if tiny
    minlen=1e-15;
else
    minlen=1e-6;
end
[xx,yy]=regexp(unrooted_newick,'(?<=:)[\d\.\-eE]*[\d\.]','split','match');
zz=str2double(yy);
if any(isnan(zz))
    yy=yy(isnan(zz));
    yy(2,:)={'...'':'};
    yy(3,:)={'''...'};
    s=[yy{[2,1,3],:}];
    error('Corrupt NEWICK string - malformed branch lengths (%s).',s);
end
tf=abs(zz)<minlen;
if any(tf) && ~tiny
    s=sprintf('Warning: %u tiny branch lengths (<10^-6) were converted to 0. (Override with ''t''.)',sum(tf));
    warning(s);
    mad_stats.messages{end+1}=s;
end
zz(tf)=0;
if any(zz<0)
    error('Cowardly refusing to root trees with negative branch lengths.')
end
if sum(zz>0)<3
    error('Cowardly refusing to root trees with less than 3 positive branch lengths.')
end
yy(tf)={'0.0'};
xx(2,1:end-1)=yy;
unrooted_newick=[xx{:}];
unrooted_newick=resolvepoly(unrooted_newick);
t=str2tree(unrooted_newick);
% end
if t.notu<3
    error('Cowardly refusing to root trees with less than 3 OTUs.');
end

%-----
%init
err='Oops..., please report this error to giddy.landan@gmail.com';
if ~all(isfinite(t.len))
   error('Corrupt NEWICK string - undefined or infinit branch lengths.');
end
if any(t.len<0)
    error('Cowardly refusing to root trees with negtive branch lengths.');
end
dis=tree2dis(t,'I'); % Distances between nodes, including internal nodes
%tip polytomies
keep=zeros(t.notu,1);
rr={};
for i=1:t.notu
    if keep(i)>0
        continue
    end
    ii=find(dis.dis(i,1:t.notu)==0);
    if numel(ii)==1
        continue
    end
    keep(ii(2:end))=i;
    dis.dis(ii(2:end),:)=nan;
    dis.dis(:,ii(2:end))=nan;
end
ntipp=sum(keep>0);
enotu=sum(keep==0);

if ntipp>0
    s=sprintf('Squeezing tip polytomies (%u OTUs, %u redundant tips, %u effective OTUs).',...
        t.notu,ntipp,enotu);
    warning(s);
    mad_stats.messages{end+1}=s;
    if enotu<3
        error('Cowardly refusing to root trees with less than 3 effective OTUs.');
    end
end
dev2=nan(t.n_branches,2,2); % Store deviations associated to each branch, leading to a local root
np=nchoosek(enotu,2); % Number of otu pairs
%-----
%---- START part one: Loop through internal nodes to compute the deviations associated to them, as local roots
for inode=t.notu+1:t.full.n_nodes(1)
    bid=t.full.node2brn(inode,:); % Get the 3 branches connected to 'inode'
    sp3=t.full.split(bid,:); % Get the 3 splits corresponding to 'bid'
    sp3(sp3(:,inode),:)=~sp3(sp3(:,inode),:); % Standardize the coding of nodes, so the cluster of the distal otus are coded as 1 and the cluster of proximal otus and 'inode' itself are coded as 0
    % Distal otus are those in the other side of the split, relative to 'inode'
    sp3=sp3(:,1:t.notu); % Remove internal nodes from splits
    d=dis.dis(1:t.notu,inode); % Distances from 'inode' to distal otus
    stt2=zeros(3,2); % Store the deviations for each pair of split
    for i=1:3 % Form pairs with the 3 splits
        da=d(sp3(2,:)); % Distances from distal otus of split 2 relative to 'inode'
        db=d(sp3(3,:)); % Distances from distal otus of split 3 relative to 'inode'
        dda=repmat(da,1,numel(db));
        ddb=repmat(db,1,numel(da))';
        dev=2.*dda./(dda+ddb) -1 ;% Relative deviations
        dev(isnan(dev))=[];
        stt2(i,1:2)=[numel(dev),sum(dev(:).^2)]; % Number of pairs, sum of squared relative deviations
        sp3=sp3([2,3,1],:); % Iterate the splits pair
    end
    nid=t.full.brn2node(bid,:); % Get nodes connected to 'bid'
    side=1+(nid(:,2)==inode); % Get the direction of the branches relative to the local root 'inode', as stored in the tree structure
    for i=1:3
        dev2(bid(i),1:2,side(i))=stt2(i,:); % The deviations implied when bid(i) leads to the local root 'inode'
    end
end
%------------------------------- END part one -----------------------------

%-------- START part two: Collect the deviations associated to each global
% root. This part does not include root expectation. It comes later ---------------
locrootnode=nan(t.n_branches,2); %%% store the position of the estimated root node on each branch of the tree, relative to 'inode'
cv=nan(t.n_branches,1);
Rb=nan(t.n_branches,1);
devt=nan(t.n_branches,2);
ronnode=false(t.full.n_nodes(1));
for bid=1:t.n_branches
    if t.len(bid)==0
        continue
    end
    nid=t.full.brn2node(bid,:);
    inode=nid(1);
    sp=t.full.split(bid,:);
    sp(sp(:,inode),:)=~sp(sp(:,inode),:);
    %one-sided pairs
    bdis=dis.brn(:,inode);
    bdis2root=bdis(t.full.brn2node);
    brn_anc_first_node=bdis2root(:,2)<bdis2root(:,1);
    brn_anc_first_node(bid)=true;
    brn_on_first_side=~sp(t.full.brn2node(:,1))'; % boolean of the branches on the same side of 'inode'
    d2x2{1,1}=dev2(brn_on_first_side & brn_anc_first_node,:,1); % Get the branches on the same side of 'inode' and uncorrectly directed in the tree structure (relative to the global root)
    d2x2{1,2}=dev2((~brn_on_first_side) & brn_anc_first_node,:,1);% Get the branches on the other side of 'inode' and uncorrectly directed in the tree structure
    brn_on_first_side(bid)=false;
    brn_anc_first_node(bid)=false;
    d2x2{2,1}=dev2(brn_on_first_side & ~brn_anc_first_node,:,2); % Get the branches on the same side of 'inode' and correctly directed in the tree structure
    d2x2{2,2}=dev2((~brn_on_first_side) & ~brn_anc_first_node,:,2); % Get the branches on the other side of 'inode' and correctly directed in the tree structure
    ancd2=cat(1,d2x2{:});
    n=size(ancd2,1);
    if n~=t.n_branches+1
        error(err);
    end
    %rooting: traversing pairs
    sp=sp(1:t.notu);
    mlen=t.len(bid);
    d=dis.dis(1:t.notu,inode);
    d1=d(~sp);
    d2=d(sp);
    dd1=repmat(d1,1,numel(d2));
    dd2=repmat(d2,1,numel(d1))';
    %-------------
    d12=(dd1+dd2);
    d2obc=d12(:).^-2;
    no=(d12(:)-2.*dd1(:)).*d2obc;
    a=nansum(no)./2./mlen./nansum(d2obc) ;% Relative deviations
    if isnan(a)
        dbstack
        keyboard
        continue
    end
    a=max(0,min(a,1));
    loc2inode=a*mlen;
    if abs(loc2inode)<minlen
        loc2inode=0;
    elseif abs(mlen-loc2inode)<minlen
        loc2inode=mlen;
    end
    
    if loc2inode==0
        if ronnode(nid(1))
            Rb(bid)=nan;
            continue
        end
        ronnode(nid(1))=true;
    end
    if loc2inode==mlen
        if ronnode(nid(2))
            Rb(bid)=nan;
            continue
        end
        ronnode(nid(2))=true;
    end

    %------------------
    locrootnode(bid,1:2)=[loc2inode,mlen-loc2inode];
    z=[d1+loc2inode;d2-loc2inode];
    cv(bid)=nanstd(z)/nanmean(z);
    dev=2.*(dd1+loc2inode)./d12 -1 ;% Relative deviations
    mad_info.type='Relative Deviations';
    dev(isnan(dev))=[];
    ancd2(end+1,:)=[numel(dev),sum(dev(:).^2)];
    devt(bid,:)=[numel(dev),sum(dev(:).^2)];
    ancd2(isnan(ancd2(:,1)),:)=[];
    if size(ancd2,1)~=t.notu-1
        error(err);
    end
    d2=sum(ancd2);
    if d2(1)~=np
        error(err);
    end
    Rb(bid)=sqrt(d2(2)/d2(1));%global mean
end
%------ END part two ---------------------------------------------------
% outputs:
madtol=1.0001;
mad_info.unrooted_newick=unrooted_newick;
mad_info.unrooted_tree=t;
mad_info.d_root_ij=locrootnode;
mad_info.r=Rb;
sRb=sort(Rb);
m=sRb(1);
ii=find(Rb<=m.*madtol);
mad_stats.nRoots=numel(ii);
if numel(mad_stats.nRoots)==1
    mad_info.ai=sRb(1)./sRb(2);
else
    mad_info.ai=1.0;
end
mad_info.cv=cv.*100;
mad_stats.root_AD=m;
mad_stats.root_AI=mad_info.ai;
mad_stats.clock_CV=mad_info.cv(ii);
mad_info.nbest=numel(ii);
mad_info.best_branch=ii;
mad_info.dev2=dev2;
mad_info.devt=devt;
mad_info.dist=dis.dis;

if numel(ii)>1
    s=sprintf('Tied root positions, %u NEWICK strings.',numel(ii));
    mad_stats.messages{end+1}=s;
    warning(s);
end
nwk1=repmat({''},1,mad_stats.nRoots);
%newick string including roots and support values
for i=1:mad_stats.nRoots
    rbid=ii(i);
    inode=t.full.brn2node(rbid,1);
    sp=t.full.split(rbid,:);
    sp(sp(:,inode),:)=~sp(sp(:,inode),:);
    sp=~sp(1:t.notu);
    rblen=mad_info.d_root_ij(rbid,:);
    if any(rblen==0)
        s='Root is polytomous';
        mad_stats.messages{end+1}=s;
        warning(s);
    end
    rt=tree2root(t,sp);
    dis=tree2dis(rt,'IR');
%     z=sum(rblen)-sum(rt.root.full.join_len(end,:));
%     if abs(z)>eps*100
%         error(err);
%     end
    bb=rt.root.full.node2brn(end,1:2);
    rsp=rt.root.full.split(bb,:);
    rsp(rsp(:,end),:)=~rsp(rsp(:,end),:);
    rsp=rsp(:,1:t.notu);
    [~,ij]=ismember(sp,rsp,'rows');
    if ij==0
        error(err);
    end
    if ij==2
        rblen=fliplr(rblen);
    end
    rt.root.full.join_len(end,:)=rblen;
    nwk1{1,i}=tree2str(rt,'R');
    s=sprintf('[MAD=%#5.3f_AI=%#5.3f_CCV=%#.3g%%_N=%u/%u]',...
        mad_stats.root_AD,mad_stats.root_AI,mad_info.cv(rbid),i,mad_stats.nRoots);
    mad_stats.messages{end+1}=s;
    if astt
        nwk1{2,i}=s;
    else
        nwk1{2,i}='';
    end
    
end
% if isempty(rr)
    rooted_newick=sprintf('%s%s;\n',nwk1{:});
% else
%     rooted_newick=regexprep(nwk1,rr(:,2),rr(:,3));
% end
mad_stats.messages{end+1}='- Please cite DOI:10.1038/s41559-017-0193';
end % function MAD
%=============================================
%=============================================
function nwstr=tree2str(tree,flags)
%----------------
ids=tree.ids;
nids=length(ids);
flags=upper(flags);
root=any(flags=='R');
%----------------
nwstr='';
%-------------
%get the right join structure
dnd=tree.root.full;
%-------------
%prepare string elements:
%node names:
nodes=cell(1,dnd.n_nodes(1)-1);
nodes(1:tree.notu)=tree.ids;
%---------
%branch lengths:
lens=cell(dnd.n_join,2);
if ~any(flags=='L')
    lens=reshape(regexp(sprintf(':%-6.5g ',dnd.join_len'),'[^ ]+','match'),2,dnd.n_join)';
    [lens{isnan(dnd.join_len)}]=deal('');%replace nans
end
%--------
%bootstrap values
bsps=cell(dnd.n_join,2);
%loop over joins:
inode=tree.notu;
for i=1:dnd.n_join-1;
    inode=inode+1;
    nodes{inode}=['(' nodes{dnd.join_nodes(i,1)} bsps{i,1} lens{i,1} ',' nodes{dnd.join_nodes(i,2)} bsps{i,2} lens{i,2} ')'];
end
%last cycle:
i=dnd.n_join;
nwstr=['(' nodes{dnd.join_nodes(i,1)} bsps{i,1} lens{i,1} ',' nodes{dnd.join_nodes(i,2)} bsps{i,2} lens{i,2} ')'];
end%function tree2str
%-----------
%=============================================
%=============================================
function tree=str2tree(ssi)
%preallocate structure array:
tree=struct('notu',[],'ids',[],'rooted',[],'root',[],'full',[],'n_branches',[],'split',[],'lengthed',[],'len',[]);
    %standardize newick strings to avoid too many \( in regexp:
    ss=regexprep(ssi,{'\s','\(','\)'},{'','<','>'});
    %-----------
    %get OTU names:
    [names,rest]=regexp(ss,'(?<=[<,])[^<>,:]*(?=[:>,])','match','split');
    n_leaf_nodes=numel(names);
    if numel(unique(names))~=n_leaf_nodes
        error('Duplicate OTU names in tree.')
    end
    tree.notu=n_leaf_nodes;
    nums=regexp(sprintf('%-u ',1:n_leaf_nodes*2-1),'\w*','match');%prepare node id strings
    %-------------
    %if input ids, reorder by input:
    otuid=1:n_leaf_nodes;
    tree.ids=names;
    
    %recode names to otuid:
    rest(2,1:n_leaf_nodes)=regexp(sprintf('<%-d>',otuid),'.*?>','match');
    ss=[rest{:}];
    %fill in nans for missing len and bsp values
    ss=regexprep(ss,'>([^:>,]*)(?=[>,])','>$1:nan');
    ss=regexprep(ss,'>:','>nan:');
    %check if rooted or unrooted
    n_nodes=sum(ss=='<');
    n_internal_nodes=n_nodes-n_leaf_nodes;
    if n_internal_nodes==n_leaf_nodes-1
        tree.rooted=1;
    elseif n_internal_nodes==n_leaf_nodes-2
        tree.rooted=0;
    else
        error('Corrupt NEWICK string ');
    end
    %-----------
    %start as if rooted:
    n_branches=2*n_leaf_nodes-2;
    n_nodes=n_branches+1;
    brn2node=nan(n_branches+1,2);
    node2brn=zeros(n_nodes,3);
    len=zeros(n_branches,1);
    
    %-----------
    %initialize loop variables:
    ibrn=[0,0];
    inode=n_leaf_nodes;
    %get pairs of nodes/branches combined by internal nodes:
    %for each match, 6 tokens: (nodeid,bsp,len)(nodeid,bsp,len)
    [match,rest]=regexp(ss,'<(\d+)>([^:]+):([^>,]+),<(\d+)>([^:]+):([^>,]+)','match','split');
    n_pairs=length(match);
    %-----------
    %loop over internal nodes
    last=0;
    while n_pairs>0
        %replace node pair with new nodeid:
        rest(2,1:n_pairs)=nums(inode+(1:n_pairs));
        ss=[rest{:}];
        v6=sscanf([match{:}],'%*c%f%*c%f%*c%f%*2c%f%*c%f%*c%f');
        if numel(v6)~=6*n_pairs
            error('corrupt NEWICK string.');
        end
        v=reshape(v6,6,n_pairs);
        for i=1:n_pairs
            inode=inode+1;
            ibrn=ibrn(2)+(1:2);
            %record branch values
            node2brn(inode,1:2)=ibrn;%two child branches
            node2brn(v([1,4],i),3)=ibrn;%add the parent branch to the child nodes
            brn2node(ibrn,1)=inode;
            brn2node(ibrn,2)=v([1,4],i);
            len(ibrn)=v([3,6],i);
        end%for i
        %get the next set of pairs:
        [match,rest]=regexp(ss,'<(\d+)>([^:]+):([^>,]+),<(\d+)>([^:]+):([^>,]+)','match','split');
        n_pairs=length(match)-last;%will evaluate to -1 on the last cycle
        if n_pairs==0 && ~tree.rooted
            last=1;%will prevent a double entry into this if block
            % unrooted trees have a trichotomy at the outer grouping
            % this translates to a remaining string with
            % one lengthed node and one unlengthed node, which
            % can appear in two oreders: <<a>:bsp:len,b> or <b,<a>:bsp:len>
            token=regexp(ss,'<(\d+)>([^:]+):([^>,]+),(\d+)>','tokens');
            if isempty(token)
                token=regexp(ss,'(\d+),<(\d+)>([^:]+):([^>,]+)>','tokens');
                token{1}=token{1}([2:4,1]);
            end
            token{1}(5:6)={'nan','0'};
            token{1}{3}(end+1)=' ';
            match{1}=sprintf(' %s',token{1}{:});
            n_pairs=1;
        end
    end%while internal nodes
    %-----------
    node2brn(end)=nan;
    brn2node=sort(brn2node,2,'descend');
    node2brn=sort(node2brn,2);
    node2brn(node2brn==0)=nan;
    n_join=n_leaf_nodes-1;
    %----------
    %initilize tree coding structure:
    dnd.n_nodes=[n_nodes,n_leaf_nodes,n_nodes-n_leaf_nodes];
    dnd.n_join=n_join;
    dnd.join_nodes=nan(n_join,2);
    dnd.join_len=dnd.join_nodes;
    %----------
    %update split and join variables
    split=false(n_branches+1,dnd.n_nodes(1));
    split(node2brn(1:n_leaf_nodes,3),1:n_leaf_nodes)=eye(n_leaf_nodes);
    node2brn(end)=n_branches+1;%dummy
    for ijoin=1:n_join
        inode=ijoin+n_leaf_nodes;
        brn3=node2brn(inode,:);
        node2=brn2node(brn3(1:2),2);
        dnd.join_nodes(ijoin,:)=node2;
        dnd.join_len(ijoin,:)=len(brn3(1:2));
        split(brn3(3),:)=split(brn3(1),:) | split(brn3(2),:);
        split(brn3(3),inode)=true;
    end
    dnd.brn2node=brn2node(1:n_branches,:);
    dnd.node2brn=node2brn;
    dnd.split=split(1:n_branches,:);
    
    %-----------
    % fill tree.root structure
    tree.root=[];
    if tree.rooted
        tree.root.full=dnd;
        tree.root.n_branches=n_branches;
        tree.root.split=dnd.split(:,1:n_leaf_nodes);
        tree.root.len=len;
    end
    
    %-----------
    %merge branches and remove root node
    dnd.n_nodes=dnd.n_nodes-[1,0,1];
    dnd.node2brn(end,:)=[];
    dnd.node2brn(dnd.node2brn==n_branches)=n_branches-1;
    dnd.brn2node(n_branches-1,:)=sort(dnd.brn2node(n_branches-[1,0],2),'descend')';
    dnd.brn2node(end,:)=[];
    len(end-1)=sum(dnd.join_len(end,:));
    len(end)=[];
    dnd.join_len(end,:)=[len(end),0];
    dnd.split(:,end)=[];
    dnd.split(end,:)=[];
    
    %-----------
    % fill tree structure
    tree.full=dnd;
    tree.n_branches=n_branches-1;
    tree.split=dnd.split(:,1:n_leaf_nodes);
    tree.lengthed=all(~isnan(len));
    tree.len=len;
    %-----------
end%function str2tree
%=============================================
%=============================================
function pw=tree2dis(tree,flags)
if nargin<2
    flags='';
else
    flags=upper(flags);
end
n=numel(tree);
pw(n)=struct('dis',[],'brn',[]);
    if any(flags=='R')
        if ~tree.rooted
            error('R flag found, but tree is not rooted');
        end
        t=tree.root;
    else
        t=tree;
    end
    if any(flags=='I')
        sp=t.full.split;
    elseif any(flags=='R')
        sp=t.full.split(:,[1:t.full.n_nodes(2),end]);
    else
        sp=t.split;
    end
    [nbrn,nnode]=size(sp);
    d=zeros(nnode);
    b=zeros(nnode);
    for i=1:nbrn
        ii=sp(i,:);
        jj=~ii;
        d(ii,jj)=d(ii,jj)+t.len(i);
        b(ii,jj)=b(ii,jj)+1;
    end
    pw.dis=d+d';
    pw.brn=b+b';
end%function tree2dis
%-----------
%-----------
%-----------
%-----------
%-----------
function [tree,rootsplit,MPflag]=tree2root(tree,outgroup_otus)
%=======================================================================
n_leaf_nodes=tree.notu;
if nargin>1
    % Backward copmatibility to version 2.0: outgroup_otus may be a split
    tree.root.type='User outgoup';
    if any(outgroup_otus==0)
        if numel(outgroup_otus)~=tree.notu
            error('corrupt input: outgroup_otus/split size ~= notu');
        end
        if ~isempty(setdiff(outgroup_otus,0:1))
            error('corrupt input: outgroup_otus/split out of range 0:1');
        end
        outgroup_otus=find(outgroup_otus);
        tree.root.type='User split';
    end
    %=======================================================================
    nog=numel(outgroup_otus);
    if nog>=tree.notu
        error('corrupt input: outgroup_otus >= notu');
    end
    if ~isempty(setdiff(outgroup_otus,1:tree.notu))
        error('corrupt input: outgroup_otus out of range 1:notu');
    end
    %find root branch:
    sp=tree.split;
    ii=find(sp(:,outgroup_otus(1))==0);
    sp(ii,:)=~sp(ii,:);
    ii=find(sum(sp(:,outgroup_otus),2)==nog);
    [m,i]=min(sum(sp(ii,:),2));
    rbrn=ii(i);%index of the root branch
    MPflag=m==nog;
    if ~MPflag
        warning('Outgroup otus are not a monophyletic clade!!!');
    end
    rootsplit=sp(rbrn,:);
    len2(1:2)=tree.len(rbrn)/2;%split the branch length equally
else
    tree.root.type='Midpoint';
    MPflag=2;
    pw=tree2dis(tree);
    [mm,ii]=max(pw.dis);
    [m,j]=max(mm);
    if m<=0%prevent crash when all branch lengths are 0
        i=1;
        j=2;
        warning('no positive pairwise distances!!!');
    else
        i=ii(j);
    end
    %find branches on the path from node i to node j
    ii=find(sum(tree.split(:,[i,j]),2)==1);
    sp=tree.split(ii,:);
    sp(sp(:,i),:)=~sp(sp(:,i),:);%make sure that node i has a 0 flag
    [bb,jj]=sort(sum(sp,2));%put the branches in the right order from j to i
    cs=cumsum(tree.len(ii(jj)));
    k=find(cs>=m/2,1);
    rbrn=ii(jj(k));%index of the root branch
    rootsplit=tree.split(rbrn,:);
    len2=cs(k)-m/2;%i-side portion of the branch length
    len2(2)=tree.len(rbrn)-len2;%j-side portion of the branch length
    % make sure that len2 and is in the same order as brn2node
    % that is, node i should be on the same side as brn2node(rbrn,1)
    if sum(tree.full.split(rbrn,[i,tree.full.brn2node(rbrn,1)]))==1
        %it isn't, so flip len:
        len2=fliplr(len2);
    end
end%if user/midpoint
%-----------
%get unrooted data:
len=tree.len;
n_branches=tree.n_branches;
dnd=tree.full;
n_nodes=dnd.n_nodes(1);
brn2node=dnd.brn2node;
node2brn=dnd.node2brn;
%-----------
%prepare new node:
newbrn=n_branches+1;
newnode=n_nodes+1;
len([newbrn,rbrn])=len2;
node2=brn2node(rbrn,:);
brn2node(rbrn,1)=newnode;
brn2node(newbrn,:)=[newnode,node2(1)];
brn2node(end+1,:)=newnode;%dummy
node2brn(node2(1),node2brn(node2(1),:)==rbrn)=newbrn;
node2brn(newnode,:)=[rbrn,newbrn,nan];
node2brn(isnan(node2brn))=newbrn+1;%dummy
%-----------
%rename nodes and branches by tracing the tree
%from the root upword.
%---
%initialize permuatation and logical indices:
permnode=1:n_leaf_nodes;
permnode(newnode)=newnode;
binnode=false(1,newnode);
permbrn=zeros(1,newbrn+1);
permbrn(newbrn)=newbrn;
permbrn(newbrn+1)=newbrn+1;%dummy
binbrn=false(1,newbrn+1);
%trace tree:
binnode(newnode)=true;%start from the root
n=1;
while n>0
    binbrn(node2brn(binnode,:))=true;%branches connected to binnode
    binnode(brn2node(binbrn,:))=true;%nodes connected to binbrn
    %rename branches from largest to smallest
    binbrn=(binbrn & ~permbrn);%reove branches that are already renamed
    n=sum(binbrn);
    permbrn(binbrn)=newbrn-(1:n);
    newbrn=newbrn-n;
    %rename nodes from largest to smallest
    binnode=(binnode & ~permnode);%reove nodes that are already renamed
    n=sum(binnode);
    permnode(binnode)=newnode-(1:n);
    newnode=newnode-n;
end%while
%--------
%rename and reorder:
permbrn(end)=nan;%rename dummy;
node2brn(permnode,:)=permbrn(node2brn);%rorder nodes and rename branches
node2brn(n_leaf_nodes+1:end,:)=sort(node2brn(n_leaf_nodes+1:end,:),2);
permbrn(end)=[];%remove dummy
brn2node(end,:)=[];%remove dumy
brn2node(permbrn,:)=permnode(brn2node);%reorder branches and rename nodes
brn2node=sort(brn2node,2,'descend');
len(permbrn)=len;
%----------
%update dnd
dnd.node2brn=node2brn;
dnd.brn2node=brn2node;
dnd.n_nodes=dnd.n_nodes+[1,0,1];
%----------
%update split and join variables
n_branches=n_branches+1;
split=false(n_branches+1,dnd.n_nodes(1));
split(node2brn(1:n_leaf_nodes,3),1:n_leaf_nodes)=eye(n_leaf_nodes);
node2brn(end)=n_branches+1;%dummy
for ijoin=1:dnd.n_join
    inode=ijoin+n_leaf_nodes;
    brn3=node2brn(inode,:);
    node2=brn2node(brn3(1:2),2);
    dnd.join_nodes(ijoin,:)=node2;
    dnd.join_len(ijoin,:)=len(brn3(1:2));
    split(brn3(3),:)=split(brn3(1),:) | split(brn3(2),:);
    split(brn3(3),inode)=true;
end%for ijoin
dnd.split=split(1:n_branches,:);
%-----------
% fill tree.root structure
tree.root.n_branches=n_branches;
tree.root.full=dnd;
tree.root.split=dnd.split(:,1:n_leaf_nodes);
tree.root.len=len;
%-----------
%merge branches and remove root node
dnd.n_nodes=dnd.n_nodes-[1,0,1];
dnd.node2brn(end,:)=[];
dnd.node2brn(dnd.node2brn==n_branches)=n_branches-1;
dnd.brn2node(n_branches-1,:)=sort(dnd.brn2node(n_branches-[1,0],2),'descend')';
dnd.brn2node(end,:)=[];
len(end)=[];
len(end)=tree.len(rbrn);
dnd.join_len(end,:)=[len(end),0];
dnd.split(:,end)=[];
dnd.split(end,:)=[];
%-----------
% fill tree structure
tree.n_branches=n_branches-1;
tree.full=dnd;
tree.split=dnd.split(:,1:n_leaf_nodes);
tree.len=len;
tree.rooted=1;
end% function tree2root
%-----------
%-----------
%-----------

function s=resolvepoly(s)
i=0;
rr={'',''};
while true
    %     c=s;
    [a,b]=regexp(s,'\([^\(\)]+\)','match','split','once');
    if strcmp(a,'')
        error('Corrupt NEWICK string.')
    end
    n=sum(a==',');
    if strcmp(b{1},'') && n<3
        if ~strcmp(b{2},';')
            error('Corrupt NEWICK string.')
        end
        break
    end
    b{2,1}=a;
    if n<2
        i=i+1;
        rr(i,1:2)={num2str(i,'@#%u#@'),a};
        b(2,1)=rr(i,1);
        s=[b{:}];
    elseif n>2 || ~strcmp(b{1,2},';')
        c=regexp(a,',','split','once');
        b{2,1}=[c{1} ',(' c{2} 'nan:0)'];
    end
    s=[b{:}];
end
rr=flipud(rr);
s=regexprep([a ';'],rr(:,1),rr(:,2));
end%function resolvepoly
%=======================================================
%==================  THE END  ==========================
%=======================================================






