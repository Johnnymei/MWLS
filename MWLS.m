%% copyright@ Xiaojun Mei, Shanghai maritime university
function estimate = MWLS(raylength,observes,NumOfSensors,Q_alpha,Q_beta)
% raylength indicates travelled length with noises under stratification propagation
% observes is the anchors with noise
% NumberOfSensors is the number of anchors
% Q_alpha is the measurement noise
% Q_beta is the anchor noise
   Operator = [-ones(NumOfSensors-1,1),eye(NumOfSensors-1,NumOfSensors-1)];
   AA = [-2*observes', ones(NumOfSensors,1)];
   bb = raylength.^2-sum(observes'.^2,2);
   w1 = inv(Q_alpha);
   u1 = (AA'*w1*AA)^-1*AA'*w1*bb;
   u = u1(1:3);
   s = observes;
   s1s = sum(observes(:,1).^2,1);
   sNs = sum(observes(:,2:end).^2,1)';
   r1s = raylength(1,:).^2;
   rNs = raylength(2:end,:).^2;
   D = zeros(NumOfSensors,NumOfSensors*3);            
      for j = 1 : NumOfSensors
          D(j,3*j-2:3*j)= (u-s(:,j))';
      end
   B = 2.*diag(raylength);
   B1 = Operator*B;
   D1 = 2*Operator*D;
   h1 = rNs-r1s-sNs+s1s;
   g1 = -2.*(s(:,2:end)-s(:,1))';
   w2 = inv(B1*Q_alpha*B1'+D1*Q_beta*D1');
   u2 = inv(g1'*w2*g1)*g1'*w2*h1;
   Gamma1 = inv(g1'*w2*g1)*g1'*w2*B1;
   Gamma2 = inv(g1'*w2*g1)*g1'*w2*D1;
   B2 = B;
      for j = 1 : NumOfSensors
          D2(j,3*j-2:3*j)= 2.*(u2-s(:,j))';
      end 
   g2 = 2.*(u2-s)';
   h2 = raylength.^2-sum(observes'.^2,2)+2.*s'*u2-sum(u2.^2,1);
   tildeG2 = [g2;eye(3)];
   tildeB2 = [B2;Gamma1];
   tildeh2 = [h2;zeros(3,1)];
   tildeD2 = [D2;Gamma2];
   w3 = inv(tildeB2*Q_alpha*tildeB2'+tildeD2*Q_beta*tildeD2');
   delta_u = -inv(tildeG2'*w3*tildeG2)*tildeG2'*w3*tildeh2;
   u3 = u2-delta_u;
   estimate = u3;
end

