function [P]=robustmeasure(s,shape,l1,n,bss,lb,ub)
% shape:         Gamma distribution parameter for shape, for D(1), shape=1: uniform distribution on the simplex
% s:             The observed sample (not in share vector form)
% n:             bootstrap sample size , typically 3
% bss            number of bootstrap samples
% lb             lower bound 5 or .05 or 10 or .1
% ub             upper bound 95 or .95 or 90 or .9
% l1             the number of samples drawn for prior - \bf{S}~D(lambda,) to calculate by simulation

%Kattuman and Ma
%This program generates sampling distributions of CoV, Gini, Entropy(H), GE, HHI, and P(S more equal than s), P(S less equal than s) from draws from Pareto
%distribution, for various alpha values 
%Also generates data for the power curve for all the indices CoV, Gini, Entropy(H), GE, HHI and P(S more equal than s), P(S less equal than s)
%for the left end of the samp dist for chosen size of the test: ie. prop from the true dist to the left of the 5% (1%)cutoff from the null hypothesized  distribution
%for the right end of the samp dist for chosen size of the test: ie. prop from the true dist to the right of the 95% (99%) cutoff from the null hypothesized  dist
%Output files
%Samp dists of different measures  for different alphas
%power test proportions at the right of the right side cutoff and left of the left side cutoff  of the samp dist for different measures

    if lb>1
        lb=lb/100;
    end
    if ub>1
        ub=ub/100;
    end

    a=1;             % "a" indexes the col number
    %pb=1;            % count
    %ps=1;            % counter for the full resampled  means

    %finding P(s more equal than S), and P(s less equal than S) 

    %resample bss times from the data vector
    for i=1:bss                  % for the number of bootstrap samples
        b=s(randperm(length(s)));   % randomly permute the smaples of size l, and choose the first n elements (bootstrap smaple size)
        x(i,:)=b(1:n);                           % x is the matrix containing bss no of bootstrap samples of size n.
    end
           
    meanx=mean(x,2);
    vasx=var(x,0,2);                   %variance with n-1
    y=sort(x,2);                       % orginal bootstrap sample (renamed) and sorted for Gini calculation
    meany=mean(y,2); 

    % convert data vector resamples to share vectors
    x=sort(x,2);                          % share vectors 
    xs= cumsum(x,2);
    xsct=bsxfun(@rdivide, x',xs(:,n)');
    x=double(xsct');
    xst=bsxfun(@rdivide, xs',xs(:,n)');
    sx=double(xst');
    sx(:,n)=1;                             % ensuring final cumulative sum element is 1.

    % use simulation if the required resample size is > 3, or if shape
    % parameter > 1; otherwise, use the geometric area 
    % if the required resample size = 3, and shape
    % parameter = 1; use geometric area
    if shape==1 &&n==3
        %Area  Approach, applied to bss resamples from the data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for c=1:bss
            a1=x(c,1);
            a2=x(c,2);
            a3=x(c,3);

            d1=sqrt((a2-a3)^2*2);
            d2=sqrt((a2-a1)^2*2);
            numwlp(c,a)=((d1+2*d2)^2-3*d2^2)/2;
            if a2==a3
                numwhp(c,a)=(3*((d2+a1*2*sqrt(2))^2-d2^2)-6*a1^2)/2;
            else if 2*(2*sqrt(2)*a1+d2)>=sqrt(2)
                numwhp(c,a)=(3*((d2+a1*2*sqrt(2))^2-d2^2)-3*(2*(2*sqrt(2)*a1+d2)-sqrt(2))^2)/2;
            else
                numwhp(c,a)=3*((d2+a1*2*sqrt(2))^2-d2^2)/2;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else

        % simulation where area is not used
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate l1 random draws from D(lambda), and translate to share vectors
        scale = 1;
        z=gamrnd(shape,scale,l1,n);
        zco=sort(z,2);
        zs= cumsum(zco,2);
        zsct=bsxfun(@rdivide, zs',zs(:,n)');
        zscosc=double(zsct');
        zscosc(:,n)=1;                 % only used for prec and succ, as comparing lorenz curves of sample with DD LC


        for c=1:bss
            for m1=1:l1
                for m2=1:n
                    msx(m1,m2)=sx(c,m2);
                end 
             end
             for m1=1:l1
                 mcheck(m1,1)=n; 
             end
     
             indexwl=zscosc>=msx;
             indexwh=zscosc<=msx;
             indexwls=sum(indexwl,2);
             indexwhs=sum(indexwh,2);
             indexwlsc=indexwls==mcheck;
             indexwhsc=indexwhs==mcheck;
             numwlp(c,a)=sum(indexwlsc)/l1;           %the probability of P(S more equal than s)  clue less - towards eq line: Schur-convex
             numwhp(c,a)=sum(indexwhsc)/l1;           %the probability of P(S less equal than s)  clue less - towards max ineq : Schur-concave
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    numwlp=sort(numwlp,1);
    numwhp=sort(numwhp,1);


    %Output collected.
    P(1,1)=mean(numwlp(1:bss,1));  %Mean of bss samples P(S more equal than s)
    P(2,1)=std(numwlp(1:bss,1));   %SD of bss samples P(S more equal than s)
    P(3,1)=numwlp(round(bss*lb),1);
    P(4,1)=numwlp(round(bss*ub),1);
    P(1,2)=mean(numwhp(1:bss,1));  %Mean of bss samples P(S less equal than s)  
    P(2,2)=std(numwhp(1:bss,1));   %SD of bss samples P(S less equal than s)
    P(3,2)=numwhp(round(bss*lb),1);
    P(4,2)=numwhp(round(bss*ub),1);

    %xlswrite('SMes.xls',numwlp)   % output: bss rows
    %xlswrite('SLes.xls',numwhp)   % output: bss rows
 
end
