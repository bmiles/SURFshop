%%%% Local membrane thinning

%%%%To do before each measuremnt -- X points

Chip = 'chip2_' %XXX
Image = '2.txt' %XXX
Duration = 180

FILE = strcat(Chip, Image)

FILEPATH= 'C:\Imperial work\Billy\Ongoing work\locallised membrane thinning\Analysis\AFM\Wafer 4 & 7\wafer 1'

FUNCTIONPATH='C:\Imperial work\Billy\Ongoing work\locallised membrane thinning\Matlab'

% 1. Extract text file for topography.

cd(FILEPATH)
% 
data = load(FILE);

% 2. Plot 2-d contour plot of topography and locate the locally thinned region

Plot1=contour(data);
title('Contour plot of relative feature height')
xlabel('X Distance/ 10^-^6 m');
ylabel('Y Distance/ 10^-^6 m');

h_ylabel=get(gca,'YLabel')
set(h_ylabel,'FontSize',18)
h_xlabel=get(gca,'XLabel')
set(h_xlabel,'FontSize',18)
set(gca,'fontsize',16)

colorbar('location','southoutside')

Figurename=['Contour_Plot_' Chip '_' Image  '.fig']
saveas(gcf,Figurename)
clear figure
close(gcf)

%3. Select region of interest!

Region=data(1:1024,1:1024); %%%XXX Change co-ordinates accordingly! Note that they are y,x! 

[C,h] =contour(Region,20)
title('Contour plot of relative feature height')
xlabel('X Distance/ 10^-^6 m');
ylabel('Y Distance/ 10^-^6 m');

h_ylabel=get(gca,'YLabel')
set(h_ylabel,'FontSize',18)
h_xlabel=get(gca,'XLabel')
set(h_xlabel,'FontSize',18)
set(gca,'fontsize',16)

colorbar('location','southoutside')

Figurename=['Local_Contour_Plot_' Chip '_' Image  '.fig']
saveas(gcf,Figurename)
clear figure
close(gcf)

% 3. Construct histogram of z value for every element within the matrix.

figure
hist(Region(:),1000)

ylabel('Count')
xlabel('Feature height/ 10^-^6 m')

Figurename=['LT_Histogram_' Chip '_' Image  '_.fig'];
saveas(gcf,Figurename)
clear figure
close(gcf)

[nelements,xcenters]=hist(Region(:),1000);

% 4. Fit 2 gaussians to histogram.

cd(FUNCTIONPATH)

X=xcenters;
Y=nelements;

x0 = [2000,0.005,0.001,12000,0.028,0.001;];

lb = [0 0 0 0 0 0];
 
options = optimset('MaxIter',10000,'TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',10000);

[x]=lsqnonlin(@Gauss2,x0,lb,[],options,X,Y);

cd(FILEPATH)

% 5. Plot 2 gaussians on chart

Plot5= scatter(xcenters,nelements)
hold on

X_Gaussian = [0:0.0001:0.1]

gaussfit = x(1).*exp(-((X-x(2)).^2).*(0.5*(x(3).^-2))) + x(4).*exp(-((X-x(5)).^2).*(0.5*(x(6).^-2)))% + x(7).*exp(-((X-x(8)).^2).*(0.5*(x(9).^-2)))%+ x(10).*exp(-((X-x(11)).^2).*(0.5*(x(12).^-2)))

p=plot(X_Gaussian, x(1).*exp(-((X_Gaussian-x(2)).*(x(3)^-1)).^2))
set(p,'Color','red','LineWidth',2)
hold on

p2=plot(X_Gaussian,x(4).*exp(-((X_Gaussian-x(5)).*(x(6)^-1)).^2))
set(p2,'Color','green','LineWidth',2)
hold on

ylabel('Count')
xlabel('Feature height/10^-^6 m')

h_ylabel=get(gca,'YLabel')
set(h_ylabel,'FontSize',18)
h_xlabel=get(gca,'XLabel')
set(h_xlabel,'FontSize',18)
set(gca,'fontsize',16)

Figurename=['LT_Histogram + Gaussian fit' Chip '_' Image  '.fig']
saveas(gcf,Figurename)
clear figure

% Extract parameters
% SD = standard deviation
% Mean = centre of peaks in micrometers

G1_count = x(1)
G1_mean = x(2)
G1_SD = x(3)
G1 = [G1_count;G1_mean;G1_SD]

G2_count = x(4)
G2_mean = x(5)
G2_SD = x(6)
G2 = [G2_count;G2_mean;G2_SD]

%%%% Etch height in nm!

Depth1 = (G2_mean-G1_mean).*10^3
Depth1_error =(((G2_SD^2)+(G1_SD^2))^0.5).*10^3

%%%% Depth 1 Etch rate

%ERate1 = Depth1/Duration
%ERate1_error = ((1/Duration)^2)*Depth1_error

%%%Calculate depth of etch using most probable values instead.

[pks,locs]=findpeaks(Y,'MINPEAKHEIGHT',150,'MINPEAKDISTANCE',200)

Peak1=[X(locs(1)),pks(1)] 
Peak2=[X(locs(2)),pks(2)]

Depth2 = (X(locs(2))-X(locs(1))).*1000

%%%% Depth 2 Etch rate

%ERate2 = Depth2/Duration

%%% Compile data

%Table1= [Duration Depth1 Depth1_error ERate1 ERate1_error Depth2 ERate2]


output_file_name = [Chip '_' Image '.mat'];
save(output_file_name)

clear all
clear figure

