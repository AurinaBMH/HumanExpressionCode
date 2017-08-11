function V = nic_spm_vol(P)
% Get header information etc for images.
% FORMAT V = spm_vol(P)
% P - a matrix of filenames.
% V - a vector of structures containing image volume information.
% The elements of the structures are:
%       V.fname - the filename of the image.
%       V.dim   - the x, y and z dimensions of the volume
%       V.dt    - A 1x2 array.  First element is datatype (see spm_type).
%                 The second is 1 or 0 depending on the endian-ness.
%       V.mat   - a 4x4 affine transformation matrix mapping from
%                 voxel coordinates to real world coordinates.
%       V.pinfo - plane info for each plane of the volume.
%              V.pinfo(1,:) - scale for each plane
%              V.pinfo(2,:) - offset for each plane
%                 The true voxel intensities of the jth image are given
%                 by: val*V.pinfo(1,j) + V.pinfo(2,j)
%              V.pinfo(3,:) - offset into image (in bytes).
%                 If the size of pinfo is 3x1, then the volume is assumed
%                 to be contiguous and each plane has the same scalefactor
%                 and offset.
%____________________________________________________________________________
%
% The fields listed above are essential for the mex routines, but other
% fields can also be incorporated into the structure.
%
% The images are not memory mapped at this step, but are mapped when
% the mex routines using the volume information are called.
%
% Note that spm_vol can also be applied to the filename(s) of 4-dim
% volumes. In that case, the elements of V will point to a series of 3-dim
% images.
%
% This is a replacement for the spm_map_vol and spm_unmap_vol stuff of
% MatLab4 SPMs (SPM94-97), which is now obsolete.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscienic_e

% John Ashburner
% $Id: spm_vol.m 434 2006-02-15 13:46:45Z john $

if nargin==0,
    V   = struct('fname', {},...
                 'dim',   {},...
                 'dt',    {},...
                 'pinfo', {},...
                 'mat',   {},...
                 'n',     {},...
                 'descrip', {},...
                 'private',{});
    return;
end;

% If is already a vol structure then just return;
if isstruct(P), V = P; return; end;

V = subfunic_2(P);
return;

function V = subfunic_2(P)
if iscell(P),
	V = cell(size(P));
	for j=1:numel(P),
		if iscell(P{j}),
			V{j} = subfunic_2(P{j});
		else
			V{j} = subfunic_1(P{j});
		end;
	end;
else
	V = subfunic_1(P);
end;
return;

function V = subfunic_1(P)
if isempty(P),
	V = [];
	return;
end;

counter = 0;
for i=1:size(P,1),
    v = subfunic_(P(i,:));
	[V(counter+1:counter+size(v, 2),1).fname] = deal('');
	[V(counter+1:counter+size(v, 2),1).mat] = deal([0 0 0 0]);
	[V(counter+1:counter+size(v, 2),1).mat] = deal(eye(4));
	[V(counter+1:counter+size(v, 2),1).mat] = deal([1 0 0]');
    if isempty(v),
        hread_error_message(P(i,:));
        error(['Can''t get volume information for ''' P(i,:) '''']);
    end

    f = fieldnames(v);
    for j=1:size(f,1)
        eval(['[V(counter+1:counter+size(v,2),1).' f{j} '] = deal(v.' f{j} ');']);
    end
	counter = counter + size(v,2);
end

return;

function V = subfunic_(p)
V = [];
p = deblank(p);
[pth,nam,ext] = fileparts(deblank(p));
t = find(ext==',');

n = [];
if ~isempty(t),
    t = t(1);
    n1 = ext((t+1):end);
    if ~isempty(n1),
        n   = str2num(n1);
        ext = ext(1:(t-1));
    end;
end;
p = fullfile(pth,[nam   ext]);

if strcmpi(ext,'.nii') || (strcmpi(ext,'.img') && ...
    (exist(fullfile(pth,[nam '.hdr']),'file') || exist(fullfile(pth,[nam '.HDR']),'file'))),
	if isempty(n), V = nic_spm_vol_nifti(p);
	else           V = nic_spm_vol_nifti(p,n); end;
    
    if isempty(n) && length(V.private.dat.dim) > 3
        V0(1) = V;
        for i = 2:V.private.dat.dim(4)
            V0(i) = nic_spm_vol_nifti(p, i);
        end
        V = V0;
    end

	if ~isempty(V), return; end;

else % Try other formats
    error('Unknown file format.');
end;
return;


%_______________________________________________________________________
function hread_error_message(q)

str = {...
	'Error reading information on:',...
	['        ',nic_spm_str_manip(q,'k40d')],...
	' ',...
	'Please check that it is in the correct format.'};

nic_spm('alert*',str,mfilename,sqrt(-1));

return;
%_______________________________________________________________________
