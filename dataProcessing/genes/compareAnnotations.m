% check if entrezID match probe ID
function [EntrezIDnomatch,ProbeIDnomatch] = compareAnnotations(ProbeName, EntrezID, GeneID, GeneSymbol, GeneName, ProbeID)

load('annotations20150612.mat')
load('initialData.mat')

% remove CUST probes - there will be no cust probes
% Remove all CUST probes (assign NaN values for all custom probes)
% fprintf(1,'Removing CUST probes\n')
% cust = strfind(ProbeName, 'CUST');
% remInd = find(~cellfun(@isempty,cust));
% fprintf(1,'%d CUST probes removed\n', length(remInd))
% ProbeName(remInd) = {NaN};
% ProbeID(remInd) = NaN;
% %fprintf(1,'Removing CUST probes\n')
% %fprintf(1,'Removing %d CUST probes\n', sum(isnan(ProbeID)))
% 
% ProbeName(isnan(ProbeID)) = [];
% EntrezID(isnan(ProbeID)) = [];
% GeneID(isnan(ProbeID)) = [];
% GeneSymbol(isnan(ProbeID)) = [];
% GeneName(isnan(ProbeID)) = [];
% ProbeID(isnan(ProbeID)) = []; 

%doesMatch = zeros(length(ProbeName), 1); 

[probes, iallen, iannot] = intersect(ProbeName, annotations20150612.ProbeID);
%getEntrez allen
entrezAllen = EntrezID(iallen); 
entrezAnnot = annotations20150612.EntrezGeneID(iannot); 

doesMatch = entrezAllen==entrezAnnot; 
EntrezIDnomatchAllen = entrezAllen(doesMatch==0); 
ProbeIDnomatchAllen = probes{doesMatch==0};
EntrezIDnomatchAnnot = entrezAnnot(doesMatch==0);


% for probe=1:length(ProbeName)
%     
%     nameAllen = ProbeName{probe};
%     entrezAllen = EntrezID(probe);
%     % for each probe in allen data find a probe in the database and check
%     % if corresponding entrezIDs match
%     out = cellfun(@(x)strcmp(nameAllen, x),annotations20150612.ProbeID,'un',0);
%     outMat = cell2mat(out);
%     ind = find(outMat==1);
%     if isempty(ind)
%         doesMatch(probe) = NaN;
%     else
%         
%         % take entrezID from there and compare to allenEntrezID
%         entrezAgilent = annotations20150612.EntrezGeneID(ind);
%         if isnan(entrezAgilent)
%             
%             doesMatch(probe) = NaN;
%         else
%             
%             doesMatch(probe) = double(entrezAllen==entrezAgilent);
%             
%         end
%     end
%    
% end
% keep those that don't match
% EntrezIDnomatch = EntrezID(doesMatch==0); 
% ProbeIDnomatch = ProbeName{doesMatch==0}; 

end

% manually check those that are 0.