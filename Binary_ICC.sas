/*******************************************************************************
    NAME    : Binary_ICC.sas
    TITLE   : ICC calculaiton for binary data
    PRODUCT : SAS R9.4
    AUTHOR  : Yosuke Inaba
*******************************************************************************/
%Macro Binary_ICC(Input=, Cluster=, Response=, Output=, Method=, CI=N, iter=50);

%* Count methods to execute;
%let cnt=%eval(%sysfunc(count(&Method.,%str( ))) + 1);

%* calculate each method ;
%do J=1 %to &cnt.;
%LET CI2=Y;
%let _Method_=%scan(&Method.,&J.,%str( ));
%*------------------------------------------------------------------------------;
%*     calculate ICC; 
%*------------------------------------------------------------------------------;
%if &_Method_.^=glmm %then %do;
proc iml;

start binary_ICC(cluster,response,method);
n=nrow(response);
u = unique(cluster);
k=ncol(u);
Y_i = j(1,k);
n_i = j(1,k);
do i = 1 to k;
    idx = loc(cluster=u[i]);
    Y_i[i] = sum(response[idx]);
    n_i[i] = nrow(response[idx]);
end;
n_0  = (1/(k-1))#(N - ((n_i##2)/N)[+]);
MS_b = (1/(k-1))#((Y_i##2/n_i)[+] - (1/N)#((Y_i[+])##2));
MS_w = (1/(N-k))#(Y_i[+] - ((Y_i##2)/n_i)[+]);
pi_i=Y_i/n_i;
result=j(1,1,.);
%* The Analysis of Variance Estimator;
if method="aov" then do;
	if Y_i[+]^=0 then do;
    rho_aov = (MS_b - MS_w)/(MS_b + (n_0 - 1)#MS_w);
    result = rho_aov;
	end;
end;

else if method="aovs" then do;
	if Y_i[+]^=0 then do;
    MS_b = (1/(k))#((Y_i##2/n_i)[+] - (1/N)#((Y_i[+])##2));
    rho_aov = (MS_b - MS_w)/(MS_b + (n_0 - 1)#MS_w);
    result = rho_aov;
	end;
end;

%*Moment Estimator;
else if method="keq" then do;
    w_i =j(1,k,1/k);
    pi_w = (w_i#pi_i)[+];
    S_w = (w_i#((pi_i - pi_w)##2))[+];
    rho_keq=(S_w - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
    result = rho_keq;
end;

else if method="kpr" then do;
    w_i =n_i/N;
    pi_w = (w_i#pi_i)[+];
    S_w = (w_i#((pi_i - pi_w)##2))[+];
    rho_keq=(S_w - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
    result = rho_keq;
end;

else if method="keqs" then do;
    w_i =j(1,k,1/k);
    pi_w = (w_i#pi_i)[+];
    S_w = (w_i#((pi_i - pi_w)##2))[+];
    S_wn = (k - 1)#S_w/k;
    rho_keq=(S_wn - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
    result = rho_keq;
end;

else if method="kprs" then do;
    w_i =n_i/N;
    pi_w = (w_i#pi_i)[+];
    S_w = (w_i#((pi_i - pi_w)##2))[+];
    S_wn = (k - 1)#S_w/k;
    rho_keq=(S_wn - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
    result = rho_keq;
end;

else if method="w" then do;
    do j=1 to 10000;
    if j=1 then rho_w=0.02;
    else do;
        _w_i_sum=(n_i/(1 + rho_w#(n_i-1)))[+];
        w_i =(n_i/(1 + rho_w#(n_i-1)))/_w_i_sum;
         pi_w = (w_i#pi_i)[+];
        S_w = (w_i#((pi_i - pi_w)##2))[+];
        rho_w=(S_w - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
    end;
    end;
    result = rho_w;
end;

else if method="ws" then do;
    do j=1 to 10000;
        if j=1 then 
            rho_ws=0.02;
        else do;
            _w_i_sum=(n_i/(1 + rho_ws#(n_i-1)))[+];
            w_i =(n_i/(1 + rho_ws#(n_i-1)))/_w_i_sum;
            pi_w = (w_i#pi_i)[+];
            S_w = (w_i#((pi_i - pi_w)##2))[+];
            S_wn = (k - 1)#S_w/k;
            rho_ws=(S_wn - (pi_w#(1-pi_w)#((w_i#(1 - w_i))/n_i)[+]))/(pi_w#(1 - pi_w)#((w_i#(1-w_i))[+]) - (w_i#(1 - w_i)/n_i)[+]);
        end;
    end;
    result = rho_ws;
end;

else if method="stab" then do;
	if Y_i[+]^=0 then do;
    kappa = 0.45;
    p = y_i[+]/n_i[+];
    w_i = n_i/N;
    pi_i = Y_i/n_i;
    pi_w = (w_i # pi_i)[+];
    s_w = (w_i # (pi_i - pi_w)##2)[+];
    rho_stab = (1/(n_0 - 1)) # ((N # s_w)/((k - 1) # p # (1 -  p)) + kappa - 1);
    result = rho_stab;
	end;
end;

else if method="ub" then do;

	if Y_i[+]^=0 then do;
	rho_ub = 1 - (N # n_0 # (k - 1) # Ms_w)/(Y_i[+] # (n_0 #
	(k - 1) - y_i[+]) + (y_i##2)[+]);
	result=rho_ub;
	end;
	end;

else if method="fc" then do;
	if Y_i[+]^=0 then do;
    pi_io=(Y_i[+])/N;
    rho_fc = 1 - (1/((N - k) # pi_io * (1 - pi_io)) # (Y_i # (n_i - Y_i)/n_i)[+]);
	result=rho_fc;
	end;
	end;

else if method="mak" then do;
	if Y_i[+]^=0 then do;
    pi_io=(Y_i[+])/N;
    rho_mak = 1 - (k - 1) # ((y_i # (n_i - Y_i))/(n_i # (n_i - 1)))[+]
        /((Y_i##2/n_i##2)[+] + (Y_i/n_i)[+] # (k - 1 - (Y_i/n_i)[+]));
	result=rho_mak;
	end;
	end;
%* Estimators Based on Diret Calculation of Correlation Within Each Group;
else if method="peq" then do;

	if Y_i[+]^=0 then do;
    mu_peq = ((n_i - 1) # Y_i)[+]/((n_i - 1) # n_i)[+];
    rho_peq = (1/(mu_peq # (1 - mu_peq))) # ((Y_i # (Y_i - 1))[+]/(n_i # (n_i - 1))[+] - mu_peq##2);
    result=rho_peq;
    end;
end;

else if method="pgp" then do;
	if Y_i[+]^=0 then do;
    mu_pgp = (Y_i/n_i)[+]/k;
    rho_pgp = (1/(mu_pgp # (1 - mu_pgp))) # (((Y_i # (Y_i - 1))/(n_i # (n_i - 1)))[+]/k - mu_pgp##2);
	result=rho_pgp;
	end;
end;

else if method="ppr" then do;
	if Y_i[+]^=0 then do;
    mu_ppr = Y_i[+]/N;
    rho_ppr = (1/(mu_ppr # (1 - mu_ppr))) # ((Y_i # (Y_i - 1)/(n_i - 1))[+]/N - mu_ppr##2);
	result=rho_ppr;
	end;
end;

return(result);

finish binary_ICC;

%* Store modules to work library;
reset storage=Binary_ICC;
store module=_all_;

quit;

proc iml;

%*Read modules;
reset storage=Binary_ICC;
load module=_all_;

%*Convert input dataset to matrix;
use &Input.;
    read all into cluster var{&Cluster.};
    read all into response var{&response.};
close &Input.;

%*Calculate ICC;
result = binary_ICC(cluster,response,"&_Method_.");

%*Output ICC value to temporary dataset;
CREATE  _mid_result1 FROM result[colname='rho_hat'];
append from result;

quit;

%* Output Warning to log if ICC is not estimable;
data _mid_result1;
    set _mid_result1;
%*    if rho_hat<-1 or 1<rho_hat then do;
%*        call missing(rho_hat);
%*        a="WARN"||"ING :  ICC is not estimable by method &_Method_..";
%*        put a; 
%*        call symputx("CI2","N");
%*    end;
%*    drop a;
run;
%end;

%else %do;
%* Calculate ICC by GLMM;

proc nlmixed data=&Input. ;
	ods table AdditionalEstimates = _mid_estimate;

  parms l_int=0.25 Sd1=0.5;
  exprhs = exp(l_int + r1);
  pred = exprhs/(1+exprhs);
  Var1=Sd1**2;
  
  model &Response. ~ binary(pred); 
  random r1 ~ normal(0,var1) subject=&Cluster.; 

  estimate 'ICC' var1/(var1 + (constant('PI')**2)/3); 

run;

data _mid_result1(keep=rho_hat);
    set _mid_estimate(rename=(estimate=rho_hat) keep=estimate);
run;

%end;

%*------------------------------------------------------------------------------;
%*     Construct 95% Confidential Interval by bootstraping; 
%*------------------------------------------------------------------------------;

%if &CI.=Y and &CI2.=Y  %then %do;
%*------  generate bootstrap sample using proc surveyselect ------;
proc surveyselect data=&Input. NOPRINT seed=1
     out=_mid_BootSSFreq(rename=(Replicate=SampleID))
     method=urs
     samprate=1
    reps=&iter.;
run;

%if &_Method_.^=glmm %then %do;
proc iml;
    %*Read modules;
    reset storage=Binary_ICC;
    load module=_all_;

    iter=&iter.;

    res=j(iter,1,0);

    do I=1 to iter;

        %*Convert dataset to matrix;
        use _mid_BootSSFreq;
            read all into cluster var{&cluster.} where(SampleID=I);
            read all into response var{&response.} where(SampleID=I);
        close _mid_BootSSFreq;

        %*Calculate ICC;
        res[I,1] = binary_ICC(cluster,response,"&_Method_.");

    end;

    %*Calculate 2.5 and 97.5 Percentiles;
    Prob = {2.5,  97.5} / 100; 
    call qntl(Pctls, res, Prob);

    %* Convert matrix to dataset;
    Pctls=Pctls`;
    varNames = {"boot_lower" "boot_upper"};
    CREATE  _mid_result2 FROM Pctls[c=varNames];
    append from Pctls;
quit;
%end;

%else %do;
%let CI2=N;
%end;

%end;
%*------  Output result ------;
data _mid_mresult_&J.;
    length method $20;
    method=symget("_Method_");
    set _mid_result1;

    %* If set CI to Y and method is not GLMM then store 95%CI;
    %if &CI.=Y and &CI2.=Y %then %do;
        if _n_=1 then set _mid_result2;
    %end;

run;

%*------  delete temporary result datasets ------;
proc datasets nolist;
    delete _mid_result:;
quit;
%end; 

%*------  Create output dataset by union result of each method ------;
data &Output.;
    set _mid_mresult_:;
run;

%*------------------------------------------------------------------------------;
%*     Delete temporary datasets; 
%*------------------------------------------------------------------------------;
proc datasets nolist;
    delete _mid_:;
quit;

%Mend Binary_ICC;
%*-- Macro end ----------------------------------------------------------------;
