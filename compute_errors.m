function res_err=compute_errors(A,R,e);
    tmp=sum(sum((A-R).^2));
    res_err=1-(tmp/e);