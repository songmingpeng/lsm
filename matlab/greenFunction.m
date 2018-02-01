% Set up arrays of X and Y coordinates
x=[1:300]; 
y=[1:300]; 
[X,Y] = meshgrid(x,y);
k=.05;
% Compute the Green¡¯s Function for k=1;
R=sqrt(X.^2+Y.^2);G=cos(k*R);
imagesc(G)
% Causal Green function if polarity=1; AcausalGreen Function if polarity=-1
polarity=1;
% Compute Movie of Propagating Wave as time increases 
for j=1:100;
    G=cos(k*R-polarity*j*.1);
    imagesc(G)
    title(['Green Function t = ',num2str(j*.01)]);
    colorbar
    xlabel('X (m)');
    ylabel('Y (m)');
    pause(.01)
end
