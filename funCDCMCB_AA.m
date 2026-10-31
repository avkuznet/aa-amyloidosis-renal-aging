function dydt = funCDCMCB_AA(t,y,hS,Akidneys,Vkidneys,Sfree,kcl,T12S, ...
k1,k2,T12A,T12B,theta12B,T12I)

% y1=[S]kidney
% y2=[A]
% y3=[B]
% y4=[I]
dydt = zeros(4,1);

dydt(1) = hS*(Akidneys/Vkidneys)*(Sfree-y(1))-kcl*y(1)-y(1)*log(2)/T12S;
dydt(2) = kcl*y(1)-k1*y(2)-k2*y(2)*y(3)-y(2)*log(2)/T12A;
dydt(3) = k1*y(2)+k2*y(2)*y(3)-y(3)*log(2)/T12B-y(3)*log(2)/theta12B;
dydt(4) = y(3)*log(2)/theta12B-y(4)*log(2)/T12I;

end