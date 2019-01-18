clear all;
clc;
% Load test data
region = 'salinasA_corrected';
class_namelist = {'class1', 'class2', 'class3', 'class4', 'class5', 'class6', 'noClass'};

% region = 'paviaU';
% class_namelist = {'Asphalt', 'BareSoil', 'Bitumen', 'Gravel', 'Meadows', 'PaintedMetalSheets', 'SelfBlockingBricks', 'Shadows', 'Trees'};

no_of_class = length(class_namelist);
for i = 1 : no_of_class
    TestData{i} = xlsread(strcat(class_namelist{i}, 'Location.xlsx'));
end
%% build confusion matrix

load('classified_image.mat');
cm = zeros(no_of_class + 1);
for i = 1 : no_of_class
    cm(i,1:no_of_class + 1) = hist(classified_image(sub2ind(size(classified_image), TestData{i}(:,1), TestData{i}(:,2))),[1:no_of_class + 1]);
end

if exist('dd'); else dd = 0; end
fpx = 1;
fpy = 1;
[nr,nc] = size(cm);
if nr ~= nc
    disp('Function Kappa: Confusion matrix must be square'); return
else
    if dd; disp(['Number of classes: ', int2str(nr)]); end
end

% marginal totals   %row sum     %Column sum
n = sum(cm(:)); cs = sum(cm); rs = sum(cm');
if dd; disp(['Number of observations: ', int2str(n)]); end
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Number of observations';
fpx = fpx + 1;
info{fpy,fpx} = int2str(n);
% naming the title of each row and column
info{fpy,7} = 'N-UA->Nive user accuracy';
info{fpy,8} = 'N-UE->Nive user error';
info{fpy,9} = 'SE-N-UA->standard error of nive user accuracy';
info{fpy,10} = 'K-UA->kappa user accuracy';
info{fpy,11} = 'SE-K-UA->standard error of kappa user accuracy';
info{fpy+1,7} = 'N-PA->Nive producer accuracy';
info{fpy+1,8} = 'N-PE->Nive producer error';
info{fpy+1,9} = 'SE-N-PA->standard error of nive producer accuracy';
info{fpy+1,10} = 'K-PA->kappa producer accuracy';
info{fpy+1,11} = 'SE-K-PA->standard error of kappa producer accuracy';
info{fpy+1,12} = 'T-PA-> Tau producer accuracy';
info{fpy+1,13} = 'SE-T-PA->standard error of Tau producer accuracy';

cm = cm(find(sum(cm)), find(sum(cm')));
[nr,nc] = size(cm);
if nr ~= nc
    disp('Function Kappa: Adjusted confusion matrix must be square'); return
else
    if dd; disp(['Number of non-zero classes: ', int2str(nr)]); end
    fpy = fpy + 1;
    fpx = 1;
    info{fpy,fpx} = 'Number of non-zero classes';
    fpx = fpx + 1;
    info{fpy,fpx} = int2str(nr);
end


% ================== Confusion Matrix Summary ================================================

n = sum(cm(:)); cs = sum(cm); rs = sum(cm');
% Confusion matrix with marginal proportions
cmx = [[cm; cs./n], [rs./n, n]'];
% overall accuracy
dsum = sum(diag(cm));
oa = dsum/n;
if dd disp(['Overall accuracy: ', num2str(oa), '; error: ', num2str(1-oa)]); end
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Overall Accuracy';
fpx = fpx + 1;
info{fpy,fpx} = num2str(oa);
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Overall Error';
fpx = fpx + 1;
info{fpy,fpx} = num2str(1 - oa);
% marginal accuracy
ua = diag(cm)'./sum(cm');%user accuracy
pa = diag(cm)'./sum(cm);%producer accuracy
if dd
    disp('User''s accuracy, errors of comission: '); disp(ua);
    disp(1-ua);
    disp(['Average user''s accuracy: ', num2str(sum(ua)/nr)]);
    disp('Producer''s reliability, errors of omission: '); disp(pa);
    disp(1-pa);
    disp(['Average producer''s reliability: ', num2str(sum(pa)/nr)]);
end
cmxx = [[cmx; [pa, mean(pa)]], [ua, mean(ua), oa]'];
if dd
    disp('Confusion matrix with marginal proportions and accuracies:');
    disp(cmxx);
end

ue = sqrt(ua.*(1 - ua)./rs ); % sigma user's accuracy
pe = sqrt(pa.*(1 - pa)./cs ); % sigma producer's accuracy

% ================== Kappa coefficient =========================================================
th1 = dsum / n; th2 = (cs*rs') / n^2;
kh = (th1-th2)/(1-th2); kh %overall kappa
th3 = sum( (cs + rs) * diag(cm) ) / n^2;
th4 = 0;
for i = 1:nr
    for j = 1:nc
        th4 = th4 + (cm(i,j) * ((cs(i) + rs(j))^2));
    end
end
th4 = th4 / n^3;
if dd; disp('th1 th2 th3 th4:'); disp([th1 th2 th3 th4]); end
th1c = 1 - th1; th2c = 1 - th2;
khv = 1/n * ...
    (     ( ( th1 * th1c ) / th2c^2 ) ...
    + ( ( 2 * th1c * ((2*th1*th2) - th3) ) / th2c^3 ) ...
    + ( ( th1c^2 * ( th4 - (4 * th2^2 ) ) ) / th2c^4 ) ...
    );% variance of overall kappa
if dd
    disp(['Overall Kappa, Variance, SE, CV : ', num2str(kh), ', ', num2str(khv), ', ', num2str(sqrt(khv)), ', ', num2str(round((sqrt(khv)/kh)*1000)/10),'%']);
end

fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Overall Kappa';
fpx = fpx + 1;
info{fpy,fpx} = num2str(kh);
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'SE';
fpx = fpx + 1;
info{fpy,fpx} = num2str(sqrt(khv));
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = '95% confidence limits for overall kappa';
fpx = fpx + 1;
info{fpy,fpx} = strcat(num2str(kh-1.96*sqrt(khv)), '-', num2str(kh+1.96*sqrt(khv)));
if dd
    disp(['95% confidence limits for kappa: ', num2str(kh-1.96*sqrt(khv)), '...', num2str(kh+1.96*sqrt(khv))]);
end

%per-class kappa, user's accuracy...
p = cm./n; uap = sum(p'); pap = sum(p); dp=diag(p)';
kpu = (dp./uap - pap)./(1 - pap);

%...and its variance
t1 = uap-dp;
t2 = (pap.*uap)-dp;
t3 = dp.*(1 - uap - pap + dp);
kpuv = ( (t1./(uap.^3 .* (1-pap).^3)) .*	((t1.*t2) + t3) )./n;

%per-class kappa, producer's reliability...
kpp = (dp./pap - uap)./(1 - uap);

%...and its variance
t1 = (pap-dp);
kppv = ( (t1./(pap.^3 .* (1-uap).^3)) .*	((t1.*t2) + t3) )./n;

if dd
    disp('Per-class kappa, Var, SE, & CV, user''s accuracy: ');
    disp(kpu); disp(kpuv); disp(sqrt(kpuv)); disp(round((sqrt(kpuv)./kpu)*1000)/10);
    disp('Per-class kappa, Var, SE, & CV, producer''s reliability: ');
    disp(kpp); disp(kppv); disp(sqrt(kppv)); disp(round((sqrt(kppv)./kpp)*1000)/10);
end


% ================== Tau coefficient =========================================================

tp = sum(cm); n = sum(cm(:));dtp = diag(cm)';
numtp = (nc.*(dtp) - tp);
denmtp = (nc-1).*tp;
tpp = numtp./denmtp;% Tau coefficient per class producer

p = cm./n; uap(1:nc) = 1/nc; pap = sum(p); dp=diag(p)';
t1 = (pap-dp);
t2 = (pap.*uap)-dp;
t3 = dp.*(1 - uap - pap + dp);

% variance of tau, producer's reliability...
tpv = ( (t1./(pap.^3 .* (1-uap).^3)) .*	((t1.*t2) + t3) )./n;
tpe = sqrt(tpv);

otc = sum(numtp)/sum(denmtp);% overall tau coefficient
% variance of overall tau coefficients
pr = 1/nc ; % number of classes
po = sum(dtp)/n ;
tvar = (1/n)*(po*(1-po)/((1 - (1/nc))^2)) ;
tse = sqrt(tvar);
tcv = round((tse/otc)*1000)/10;

fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Overall Tau';
fpx = fpx + 1;
info{fpy,fpx} = num2str(otc);
fpy = fpy + 1; fpx = 1;
info{fpy,fpx} = 'Standard error of overall tau';
fpx = fpx + 1;
info{fpy,fpx} = tse;
fpy = fpy + 2; fpx = 1;


% =============== writing all the generated statistics =============

info{fpy,fpx} = 'Confusion matrix with marginal proportions and accuracies';
fpy = fpy + 1;
for i = 1:nr
    info{fpy,i+1} = strcat('class ',num2str(i));
    info{fpy+i,1} = strcat('class ',num2str(i));
end
info{fpy,nr+2} = 'RowSum';
info{fpy,nr+3} = 'N-UA';
info{fpy,nr+4} = 'N-UE';
info{fpy,nr+5} = 'SE-N-UA';
info{fpy,nr+6} = 'K-UA';
info{fpy,nr+7} = 'SE-K-UA';

info{fpy+nr+1,1} = 'ColumnSum';
info{fpy+nr+2,1} = 'N-PA';
info{fpy+nr+3,1} = 'N-PE';
info{fpy+nr+4,1} = 'SE-N-PA';
info{fpy+nr+5,1} = 'K-PA';
info{fpy+nr+6,1} = 'SE-K-PA';
info{fpy+nr+7,1} = 'T-PA';
info{fpy+nr+8,1} = 'SE-T-PA';
fpx = 2;
fpy = fpy + 1;
for i = fpy:fpy+nr-1
    for j = 2:2+nr-1
        info{i,j} = cm(i-fpy+1,j-1);
    end
end
for i = 1:nr
    info{fpy+i-1,fpx+nr} = rs(i);
    info{fpy+i-1,fpx+nr+1} = ua(i);
    info{fpy+i-1,fpx+nr+2} = 1-ua(i);
    info{fpy+i-1,fpx+nr+3} = ue(i);
    info{fpy+i-1,fpx+nr+4} = kpu(i);
    info{fpy+i-1,fpx+nr+5} = sqrt(kpuv(i));
end
fpy = fpy + nr;
info{fpy, nr + 2} = n;
for i = 1:nr
    info{fpy,fpx+i-1} = cs(i);
    info{fpy+1,fpx+i-1} = pa(i);
    info{fpy+2,fpx+i-1} = 1-pa(i);
    info{fpy+3,fpx+i-1} = pe(i);
    info{fpy+4,fpx+i-1} = kpp(i);
    info{fpy+5,fpx+i-1} = sqrt(kppv(i));
    info{fpy+6,fpx+i-1} = tpp(i);
    info{fpy+7,fpx+i-1} = tpe(i);
end

% Finally write all data to file.
xlswrite(strcat(region,'.xlsx'),info);
