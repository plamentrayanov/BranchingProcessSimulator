function Im_matrix = Im_function(Im_mass, Im_age_struct_prob)
Im_matrix=zeros(size(Im_age_struct_prob,1), 1, size(Im_mass,2));
Im_matrix(:,1,:)=mnrnd(binornd(100,Im_mass'), Im_age_struct_prob)';
end

