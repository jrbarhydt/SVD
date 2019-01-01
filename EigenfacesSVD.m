% SVD Basis Faces: Spectrum Analysis, Projection, and Reconstruction
% Johnathon R Barhydt
% 
% Yale Face Files Located at:
% http://vision.ucsd.edu/~iskwak/ExtYaleDatabase/ExtYaleB.html
% 

clear all, close all, clc

% photo directory name/location:
% uncomment to change between cropped/uncropped sample set.
files = dir( 'CroppedYale/*/*.pgm' );
%files = dir( 'yalefaces/subject*.*' );

% image size: change multiplier to maintain aspect ratio.
% 0.25 is good for speed.
multiplier = 0.25;
im_height = 192*multiplier;
im_width = 168*multiplier;

% grayscale photo not in set to reconstruct.
photo = im2double( imresize( imread('burke.jpg'), [im_height im_width]));

% reconstruction mode truncation tolerance.
tol = 1e-3;

% load face database.
A = zeros( im_height*im_width, length( files ));
for i = 1:length( files )
    image = imread( strcat( files(i).folder,'\',files(i).name) );
    A_matrix = im2double( imresize( image, [im_height im_width]));
    A(:,i) = A_matrix(:);
end

%% perform the SVD (computationally expensive).
[u,s,v] = svd(A,'econ');

%% collect singular values, normalize, determine truncation.
sing_values = diag(s);
s_norm = sing_values/sum(sing_values);
data_size = length(s_norm);
[~,trunc] = min(abs(s_norm-tol*ones(data_size,1)));

% SVD spectrum plot (raw logarithmic and normalized).
figure(1)
ax1 = subplot(2,1,1);
semilogy(ax1, sing_values, 'r.'), xlim([0 data_size])
hold on
cutoff = plot(ax1,1:data_size,sing_values(trunc)*ones(data_size,1));
hold off
title(ax1, 'Singular Value Full Spectrum')
legend(cutoff,['Truncation Cutoff at %',num2str(100*tol)])
ylabel(ax1, {'$log[s_{jj}]$'},'Interpreter','latex')
xlabel(ax1, 'index j')
ax2 = subplot(2,1,2);
bar(ax2, s_norm(1:50))
title(ax2, 'Dominant Singular Values - Normalized')
dim = [.6 0 0 .4];
str={'Normaliztion Constant:';'A = $\sum_{j=1}^{n} s_{jj}$'};
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex')
ylabel(ax2, {'$\frac{1}{A} s_{jj}$'},'Interpreter','latex')
xlabel(ax2, 'index j')

% reconstruction percentage plot
recon_percent = zeros(data_size,1);
for i = 1:data_size
   recon_percent(i) = sum(s_norm(1:i)); 
end

figure(2), area(recon_percent), xlim([0 data_size]), ylim([0 1])
title('Reconstruction Percent')
xlabel(['Number of Basis Face Vectors Used (out of ',num2str(data_size),')'])
labels = reshape(sprintf('%5.0f%%',0:10:100),6,[]).';
set(gca,'yticklabel',labels)

% prominent eigenfaces (color-reversed for viewing pleasure)
figure(3)
for i = 1:7
    ax(i)=subplot(1,7,i);
    imshow(mat2gray( reshape( -u(:,i),[im_height im_width]) ) )
end
title(ax(4),'Top 8 Basis Faces by Rank','FontSize',13)
annotation('textarrow',[.7 .2],[.1 .1],'String','Increasing Dominance')

%% facial reconstruction using k modes (computationally semi-costly).
%  note: reconstructs ALL of the faces, displaying only chosen original_face.
original_face = 149; % pick a face from the photo directory

figure(4)
k = [1 5 10 20 30 50 100 trunc data_size];
for j=1:9
    reface = u(:,1:k(j))*s(1:k(j),1:k(j))*v(:,1:k(j)).';
    ax(j) = subplot(3,3,j);
    imshow( reshape( reface(:,original_face),[im_height im_width]))
    str = ['Rank: ' num2str(k(j))];
    xlabel(str)
end
title(ax(2),'Facial Reconstruction by Rank','FontSize',13)

%% randomly selects single face to reconstruct 
%  adds noisy version along with photo not in set.
selector = 301;%randi([1,data_size]);
image = zeros(im_height*im_width,3);
image(:,1) = A(:,selector);
image(:,2) = A(:,selector)+ 0.4*randn(im_height*im_width,1); %noisy version
image(:,3) = reshape( photo, [im_height*im_width,1]);

% basis projection coefficients
b = zeros(data_size, 3);
b(:,1)=image(:,1)'*u;
b(:,2)=image(:,2)'*u;
b(:,3)=image(:,3)'*u;

% reconstruct X iteratively from basis coeffs up to truncation
X = zeros(im_height*im_width,3);
for i = 1:3
    for k = 1:trunc
        X(:,i) = X(:,i)+b(k,i)*u(:,k);
    end
end

% show originals, spectrum 'fingerprint', and reconstruction
figure(5)
for i = 1:3
    subplot(3,3,3*i-2)
    imshow( mat2gray( reshape( image(:,i),[im_height im_width])))
    subplot(3,3,3*i-1)
    bar(b(1:trunc, i)), ylim([-15 15])
    set(gca,'yticklabel',[],'xticklabel',[])
    subplot(3,3,3*i)
    imshow( mat2gray( reshape( X(:,i),[im_height im_width])))
end