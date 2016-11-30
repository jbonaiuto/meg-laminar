function trial_tc=extract_individual_trial(t, M, U, D, nverts, ntimes)
    d1 = squeeze(D(:,:,t));
    Dtrial=M*U*d1;
    trial_tc=zeros(nverts,ntimes);
    for j=1:ntimes
        trial_tc(:,j)=sum(Dtrial(:,max([1 j-25]):min([ntimes j+25])).^2,2);
    end
    %trial_tc=movingvar(Dtrial',100)';
end

