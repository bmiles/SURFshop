% 3 gaussian fit!

function gaussfit= Gauss2(x,X,Y)

A=x(1);
B=x(2);
C=x(3);
D=x(4);
E=x(5);
F=x(6);
% G=x(7);
% H=x(8);
% I=x(9);

gaussfit = (Y-(A.*exp(-((X-B).^2).*(0.5*(C.^-2)))+ D.*exp(-((X-E).^2).*(0.5*(F.^-2)))))% + G.*exp(-((X-H).^2).*(0.5*(I.^-2)))))

end