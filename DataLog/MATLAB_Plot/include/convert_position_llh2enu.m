function posenu = convert_position_llh2enu(llh_vec1, llh_vec2)
%   Convert from two lon, lat, hgt to local-level position
%
%	posenu = convert_position_llh2enu(llh_vec1, llh_vec2)
%
%    INPUTS
%	llh_vec1 = user Llh position matrix or vector in degrees and meters
%	llh_vec2 = reference Llh position matrix or vector in degrees and
%	           meters
%
%    OUTPUTS
%	posenu = local-level position matrix or vector difference between user
%            and referece position
% 
% 
[n1 m1] = size(llh_vec1);
[n2 m2] = size(llh_vec2);

posenu = zeros(n1,m1);

if n2 == 1,
    orgxyz = llh2xyz([llh_vec2(1:2)*pi/180 llh_vec2(3)]);
    for i = 1:n1,
        posxyz = llh2xyz([llh_vec1(i,1:2)*pi/180 llh_vec1(i,3)]);
        posenu(i,:) = xyz2enu(posxyz, orgxyz);
    end
elseif n2 > 1 && n1 >= n2,
    for i = 1:n2,
        orgxyz = llh2xyz([llh_vec2(i,1:2)*pi/180 llh_vec2(i,3)]);
        posxyz = llh2xyz([llh_vec1(i,1:2)*pi/180 llh_vec1(i,3)]);
        posenu(i,:) = xyz2enu(posxyz, orgxyz);
    end
else
    for i = 1:n1,
        orgxyz = llh2xyz([llh_vec2(i,1:2)*pi/180 llh_vec2(i,3)]);
        posxyz = llh2xyz([llh_vec1(i,1:2)*pi/180 llh_vec1(i,3)]);
        posenu(i,:) = xyz2enu(posxyz, orgxyz);
    end   
end