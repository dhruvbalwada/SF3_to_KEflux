%% Rotational-divergent decomposition function 
% do the integrals to calculate the decomposition to rotational and divergent part (Lindborg 2015)
% can potentially improve this by using a proper way to calculate the
% integral instead of the nansum
%clear s2rr s2dd

function [SF2rr SF2dd] = helmholtz_decompose(dist_axis, SF2ll, SF2tt)
mid_dist_axis = 0.5*(dist_axis(1:end-1)+dist_axis(2:end));
mid_diff_du = 0.5*((SF2tt(1:end-1) - SF2ll(1:end-1))+(SF2tt(2:end)-SF2ll(2:end)));

SF2rr(1) = SF2tt(1);
SF2dd(1) = SF2ll(1);
    for i =2:length(dist_axis)
        SF2rr(i) = SF2tt(i) + nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1).*diff(dist_axis(1:i)));
        SF2dd(i) = SF2ll(i) - nansum(1./mid_dist_axis(1:i-1).*mid_diff_du(1:i-1).*diff(dist_axis(1:i)));
    end
end    