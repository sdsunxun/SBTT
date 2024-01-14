function [cout,marked_img] =SBTT(varargin)
%   CSSCORNER Find corners in an image by using the scale-space concept.
%
%       CSSCORNER works by the following step:
%       1.	Apply the Canny edge detector to the gray level image and obtain a
%       binary edge-map.
%       2.	Extract the edge contours from the edge-map, fill the gaps in the
%       contours.
%       3.	Compute curvature at a low scale for each contour to retain all
%       true corners.
%       4.	All of the curvature local maxima are considered as corner
%       candidates, then rounded corners and false corners due to boundary
%       noise and details were eliminated.
%       5.  End points of line mode curve were added as corner, if they are not
%       close to the above detected corners.
%
%       Syntax :
%%       [cout,marked_img]=corner(I,C,T_angle,sig,H,L,Endpiont,Gap_size)

%       Input :
%       I -  the input image, it could be gray, color or binary image. If I is
%           empty([]), input image can be get from a open file dialog box.
%       C -  denotes the minimum ratio of major axis to minor axis of an ellipse,
%           whose vertex could be detected as a corner by proposed detector.
%           The default value is 1.5.
%       T_angle -  denotes the maximum obtuse angle that a corner can have when
%           it is detected as a true corner, default value is 162.
%       Sig -  denotes the standard deviation of the Gaussian filter when
%           computeing curvature. The default sig is 3.
%       H,L -  high and low threshold of Canny edge detector. The default value
%           is 0.35 and 0.
%       Endpoint -  a flag to control whether add the end points of a curve
%           as corner, 1 means Yes and 0 means No. The default value is 1.
%       Gap_size -  a paremeter use to fill the gaps in the contours, the gap
%           not more than gap_size were filled in this stage. The default
%           Gap_size is 1 pixels.
%
%       Output :
%       cout -  a position pair list of detected corners in the input image.
%       marked_image -  image with detected corner marked.
%
%       Examples
%       -------
%       I = imread('alumgrns.tif');
%       cout = corner(I,[],[],[],0.2);
%
%       [cout, marked_image] = corner;
%
%       cout = corner([],1.6,155);
%
%
%   Composed by Baojiang Zhong


% clc;
% warning off;
% close all;
[I,H,L,Gap_size] = parse_inputs(varargin{:});

if size(I,3)==3
    % Transform RGB image to a Gray one.
    I=rgb2gray(I);
end

tic
% Detect edges
BW=edge(I,'canny',[L,H]);
time_for_detecting_edge = toc;

tic
% Extract curves
[curve,curve_start,curve_end,curve_mode,curve_num]=extract_curve(BW,Gap_size);
time_for_extracting_curve = toc;

tic

if size(curve{1})>0


    % Detect corners
    cout  = get_corner(I,curve,curve_start,curve_end,curve_mode,curve_num);
    time_for_detecting_corner = toc;
    
        figure; imshow(BW); hold on; plot(cout(:,2),cout(:,1),'r*');

%     img=I;
%     for i=1:size(cout,1)
%         img=mark(img,cout(i,1),cout(i,2),5);
%     end
%     marked_img=img;
    
%     figure
%     clf reset
%     imshow(marked_img);
%     title('Detected corners')
else
    cout = [];
    marked_img = [];
    cd = [];
    
end
% imwrite(marked_img,'corner.jpg');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [curve,curve_start,curve_end,curve_mode,cur_num]=extract_curve(BW,Gap_size)

%   Function to extract curves from binary edge map, if the endpoint of a
%   contour is nearly connected to another endpoint, fill the gap and continue
%   the extraction. The default gap size is 1 pixles.

[L,W]=size(BW);
BW1=zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge=zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W)=BW;
[r,c]=find(BW1==1);
cur_num=0;

while size(r,1)>0
    point=[r(1),c(1)];
    cur=point;
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [~,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[cur;point];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    % Extract edge towards another direction
    point=[r(1),c(1)];
    BW1(point(1),point(2))=0;
    [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    while size(I,1)>0
        dist=(I-Gap_size-1).^2+(J-Gap_size-1).^2;
        [~,index]=min(dist);
        point=point+[I(index),J(index)]-Gap_size-1;
        cur=[point;cur];
        BW1(point(1),point(2))=0;
        [I,J]=find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
    
    if size(cur,1)>(size(BW,1)+size(BW,2))/25
        cur_num=cur_num+1;
        curve{cur_num}=cur-Gap_size;
    end
    [r,c]=find(BW1==1);
    
end

for i=1:cur_num
    curve_start(i,:)=curve{i}(1,:);
    curve_end(i,:)=curve{i}(size(curve{i},1),:);
    if (curve_start(i,1)-curve_end(i,1))^2+...
            (curve_start(i,2)-curve_end(i,2))^2<=8
        curve_mode(i,:)='loop';
    else
        curve_mode(i,:)='line';
    end
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L)=1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cout = get_corner(I,curve,curve_start,curve_end,curve_mode,curve_num)
% Control the smoothness of the curve; The larger the value, the higher the smoothness, the higher the false positives, but the more missed detections
c_smoothing = 0.3; 
% Window curvature threshold for filtering rounded corner points; The larger the value, the stronger the filtering ability
c_roundremove = 0.12; 
sig=2;
c_curvaturetherold=0.12;

        
missed_method_state=1;
round_method_state=1;
he_method_state=0;

corner_num=0;
cout=[ ];
% m is to record the number of detected corners
m(1:curve_num) = 0;


sn = 0;
for h = 1:curve_num;
    % n is the number of points on the curve.
    n_h = size(curve{h},1);
    sn = sn + n_h;
end

[gau W] = makeGFilter(sig);



        
css_out=[];
missed_method_out=[];
round_method_out=[];
he_method_out=[];

for h = 1:curve_num;
    
    
    L = size(curve{h},1);
    
    iter = ceil(c_smoothing*L);

    
    % n is the number of points on the curve.
    
    x = curve{h}(:,1);
    y = curve{h}(:,2);
    
    
    curveLen = size(x,1);    
    [xs ys W] = smoothing(x,y,curveLen,curve_mode(h,:),gau,W); % smooth the curve with Gaussian kernel
    x=xs;
    y=ys;
    
    % A is to record the CSS trajectories of corners
    A = zeros(iter,L);
    % B is to save the curvature of corner points
    B = A;
    % curvautrure of points
    kappa(1:L) = 0;
    
    xn(1:L) = 0;
    yn(1:L) = 0;
    
    first_extremum_index=[];
    first_extremum_count=0;
    
    for k = 1 : iter
        
        % in the case of a line
        if strcmp(curve_mode(h,:),'line') == 1
            
            % computing the curvature
            for i = 2 : L-1
                delta_x = (x(i+1)-x(i-1))/2;
                delta_y = (y(i+1)-y(i-1))/2;
                delta_2_x = x(i+1)-2*x(i)+x(i-1);
                delta_2_y = y(i+1)-2*y(i)+y(i-1);
                kappa(i)=(delta_x*delta_2_y-delta_y*delta_2_x)...
                    /(delta_x^2+delta_y^2)^1.5;
            end
            
            for j = 3 : L-2
                if ((abs(kappa(j)) >  abs(kappa(j-1))) && (abs(kappa(j)) > abs(kappa(j+1))))
                    A(k,j) = sign(kappa(j));
                    B(k,j) = kappa(j);
                end
            end
            
            % the beginning of gaussian smoothing
            xn(1) = 0.75*x(1)+0.25*x(2);
            yn(1) = 0.75*y(1)+0.25*y(2);
            for i = 2 : L-1
                xn(i) = 0.25*x(i-1) + 0.5*x(i) + 0.25*x(i+1);
                yn(i) = 0.25*y(i-1) + 0.5*y(i) + 0.25*y(i+1);
            end
            xn(L) = 0.25*x(L-1)+0.75*x(L);
            yn(L) = 0.25*y(L-1)+0.75*y(L);
            
            for i = 1 : L
                x(i)=xn(i);
                y(i)=yn(i);
            end
            
        end % end if; in the case of a line
        
        % in the case of a loop
        if strcmp(curve_mode(h,:),'loop') == 1
            
%         if k==41
%             k
%         end
            
            delta_x = (x(2)-x(L))/2;
            delta_y = (y(2)-y(L))/2;
            delta_2_x = x(2)-2*x(1)+x(L);
            delta_2_y = y(2)-2*y(1)+y(L);
            kappa(1)=(delta_x*delta_2_y-delta_y*delta_2_x)/(delta_x^2+delta_y^2)^1.5;
            
            % computing the curvature
            for i = 2 : L-1
                
%                 if i==140
%                 end
                delta_x = (x(i+1)-x(i-1))/2;
                delta_y = (y(i+1)-y(i-1))/2;
                delta_2_x = x(i+1)-2*x(i)+x(i-1);
                delta_2_y = y(i+1)-2*y(i)+y(i-1);
                kappa(i)=(delta_x*delta_2_y-delta_y*delta_2_x)/(delta_x^2+delta_y^2)^1.5;
            end
            
            delta_x = (x(1)-x(L-1))/2;
            delta_y = (y(1)-y(L-1))/2;
            delta_2_x = x(1)-2*x(L)+x(L-1);
            delta_2_y = y(1)-2*y(L)+y(L-1);
            kappa(L)=(delta_x*delta_2_y-delta_y*delta_2_x)/(delta_x^2+delta_y^2)^1.5;
            
            if ((abs(kappa(1)) >  abs(kappa(L))) && (abs(kappa(1)) > abs(kappa(2)))) % 角点确定
                A(k,1) = sign(kappa(1));
                B(k,1) = kappa(1);
            end
            for j = 2 : L-1
                if ((abs(kappa(j)) >  abs(kappa(j-1))) && (abs(kappa(j)) > abs(kappa(j+1)))) % 角点确定
                    A(k,j) = sign(kappa(j));  
                    B(k,j) = kappa(j);  
                end
            end
            if ((abs(kappa(L)) >  abs(kappa(L-1))) && (abs(kappa(L)) > abs(kappa(1)))) % 角点确定
                A(k,L) = sign(kappa(L));
                B(k,L) = kappa(L);
            end
            
            % Gaussian smoothing
            xn(1) = 0.25*x(L) + 0.5*x(1) + 0.25*x(2);
            yn(1) = 0.25*y(L) + 0.5*y(1) + 0.25*y(2);
            for i = 2 : L-1
                xn(i) = 0.25*x(i-1) + 0.5*x(i) + 0.25*x(i+1);
                yn(i) = 0.25*y(i-1) + 0.5*y(i) + 0.25*y(i+1);
            end
            xn(L) = 0.25*x(L-1) + 0.5*x(L) + 0.25*x(1);
            yn(L) = 0.25*y(L-1) + 0.5*y(L) + 0.25*y(1);
            
            for i = 1 : L
                x(i)=xn(i);
                y(i)=yn(i);
            end
            
            
        end % end if; in the case of a loop

%         figure;plot(y,x);
%         set(gca,'ydir','reverse')
%         axis off
    end % end for; the computing of CSS trajectory map
    

    
    for i = 1 : L
        % to locate the coner points at the highest scale
        if A(iter,i) ~= 0
            first_extremum_count = first_extremum_count+1;
            % cp(i) means the i'th corner point
            % the m'th corner locates at the i'th point
            first_extremum_index(first_extremum_count)= i;
        end
    end
    
    
    D=zeros(iter,first_extremum_count);
    if first_extremum_count > 0
        % cp_old is for relocate the lost corners
        % because cp(i) could disappear at the lowest scale
        for i = 1:first_extremum_count
            % to trace the coners back to the lowest scale
            corner_location =first_extremum_index(i);
            for k = iter-1 : -1 : 1
                % move of the trajectory is assumed to be 1 pixel
                for j = first_extremum_index(i)-1 : first_extremum_index(i)+1
                    l = j;
                    if j < 1; l = j+L; end
                    if j > L; l = j-L; end
                    if A(k,l) ~= 0
                        corner_location = l;
                    end
                    D(k,i)=corner_location;
                    
                end
                % the location of the corners is adjusted
                first_extremum_index(i) = corner_location;
            end
        end


       %% missing corner recover begin
       
            remiss_method_extremum_index=first_extremum_index;   %m_1 = m(h);
            first_extremum_count=size(first_extremum_index,2);    
            % m_1 is the number of corners for the first run
            remiss_method_extremum_count = first_extremum_count;

            % m_1 is the number of corners for the first run

            % range to a coner to relocate the new corners
            b_rec = round(iter/10);
            % search of new coners starts at a lower sclae k_rec
            k_rec = ceil(iter/25);

            for i = 1 : first_extremum_count
                for j = first_extremum_index(i)- b_rec :first_extremum_index(i) + b_rec
                    l = j;
                    if j < 1; l = j+L; end
                    if j > L; l = j-L; end
                    % for END model, relocate from k_rec;
                    % for STAIR model, relocate from 2*k_rec;
                    if l==139
                    end
                    if (A(k_rec,l)*A(k_rec,D(k_rec,i)) == 1) || (A(2*k_rec,l)*A(k_rec,D(k_rec,i)) == -1)
                        if(abs(B(k_rec,l))>c_curvaturetherold)
                            remiss_method_extremum_count = remiss_method_extremum_count +1;
                            remiss_method_extremum_index(remiss_method_extremum_count) = l;
                        end
                    end
                end
            end

            % to trace new coners back to the lowest scale
            for i = first_extremum_count+1 : remiss_method_extremum_count
                corner_location = remiss_method_extremum_index(i);
                for k = 2*k_rec-1 : -1 : 1
                    for j = remiss_method_extremum_index(i)-1 : remiss_method_extremum_index(i)+1
                        l = j;
                        if j < 1; l = j+L; end
                        if j > L; l = j-L; end
                        if A(k,l) ~= 0
                            corner_location = l;
                        end
                    end
                    remiss_method_extremum_index(i) = corner_location;
                end
            end

            % to delete the corners which are responds more than twice
            remiss_method_extremum_index = unique(remiss_method_extremum_index);
            
            extremum_index=remiss_method_extremum_index;
       %% missing corner recover end

        
        %% round corner remove begin
            % to delete round corners
            third_extremum_index=[];
            third_extremum_curvature=[];
            for i = 1:length(extremum_index);
                stable_corner_location(1:iter) = extremum_index(i);
                C(1,i) = B(1,extremum_index(i));
                % to trace the coners back to the lowest scale
                for k = 2 : iter
                    % move of the trajectory is assumed to be 1 pixel
                    for j = stable_corner_location(k-1) -1 : stable_corner_location(k-1) +1
                        l = j;
                        if j < 1; l = j+L; end
                        if j > L; l = j-L; end
                        if A(k,l) ~= 0
                            stable_corner_location(k) = l;
                            C(k,i) = B(k,l);
                        end
                    end
                end
                % when stable_corner = 0, the stability rule is used and the
                % detected corner is loated at its stable location
                stable_corner = 0;
                if stable_corner
                    T = tabulate(stable_corner_location);
                    [F,I]=max(T(:,2));
                    if F>=iter*0.5
                        dcp{h}(i)=T(I,1);
                    end
                end
            end
             front_range=floor(iter*0.1);
             behind_range=floor(iter*0.3);

            if(h==19)
            end
            for  i = 1:length(extremum_index);
                %if(size(dcp_curvature{h},1)>20)
                    if(size(C,1)>behind_range&&front_range>0)
                        
%                         i
%                         C([front_range:behind_range],i);
                        s_curvature = sum((C([front_range:behind_range],i).^2));
                        if s_curvature >= c_roundremove
                            third_extremum_index=[third_extremum_index,extremum_index(i)];
                        end
                    else
                        C
                    end
                %end
            end
            
            extremum_index=third_extremum_index;
%             for i=1:size(extremum_index,2)   
%                 coor=curve{h}(extremum_index(i),:); 
%                 round_method_out=[round_method_out;coor(1,2),coor(1,1)];
%             end

        %% round corner remove End
            count=size(third_extremum_index,2);     %%
            for j=1:count
                corner_num=corner_num+1;
                cout(corner_num,:)=curve{h}(third_extremum_index(j),:); %输出坐标
            end
            count=size(remiss_method_extremum_index,2);     %%
            for j=1:count     
                corner_num=corner_num+1;
                cout(corner_num,:)=curve{h}(remiss_method_extremum_index(j),:); %输出坐标
            end      
        
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ang=curve_tangent(cur,center)

for i=1:2
    if i==1
        curve=cur(center:-1:1,:);
    else
        curve=cur(center:size(cur,1),:);
    end
    L=size(curve,1);
    if L>3
        if sum(curve(1,:)~=curve(L,:))~=0
            M=ceil(L/2);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M,1);
            y2=curve(M,2);
            x3=curve(L,1);
            y3=curve(L,2);
        else
            M1=ceil(L/3);
            M2=ceil(2*L/3);
            x1=curve(1,1);
            y1=curve(1,2);
            x2=curve(M1,1);
            y2=curve(M1,2);
            x3=curve(M2,1);
            y3=curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2*(-y1*x2^2+y3*x2^2-y3*y1^2-y3*x1^2-y2*y3^2+x3^2*y1+y2*y1^2-y2*x3^2-y2^2*y1+y2*x1^2+y3^2*y1+y2^2*y3)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            y0 = -1/2*(x1^2*x2-x1^2*x3+y1^2*x2-y1^2*x3+x1*x3^2-x1*x2^2-x3^2*x2-y3^2*x2+x3*y2^2+x1*y3^2-x1*y2^2+x3*x2^2)/(-y1*x2+y1*x3+y3*x2+x1*y2-x1*y3-x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_direction=angle(complex(x0-x1,y0-y1));
            adjacent_direction=angle(complex(x2-x1,y2-y1));
            tangent_direction=sign(sin(adjacent_direction-radius_direction))*pi/2+radius_direction;
        end
    
    else % very short line
        tangent_direction=angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    direction(i)=tangent_direction*180/pi;
end
ang=abs(direction(1)-direction(2));

function [G W] = makeGFilter(sig);

    GaussianDieOff = .0001; 
    pw = 1:100;

    ssq = sig*sig;
    W = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
    if isempty(W)
        W = 1;  
    end
    t = (-W:W);
    gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
    G=gau/sum(gau);

function [I,H,L,Gap_size,args] = parse_inputs(varargin)

%Default experience value;
Para=[0.35,0,1];

if nargin>2
    I=varargin{1};
    for i=2:nargin
        if size(varargin{i},1)>0
            Para(i-1)=varargin{i};
        end
    end
end

if nargin==1
    I=varargin{1};
end

if nargin==2
    I=varargin{1};
    args=varargin{2};
end

if nargin==0 || size(I,1)==0
    [fname,dire]=uigetfile('*.bmp;*.jpg;*.gif','Open the image to be detected');
    I=imread([dire,fname]);
end

H=Para(1);
L=Para(2);
Gap_size=Para(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function img1=mark(img,x,y,w)

[M,N,C]=size(img);
img1=img;

if isa(img,'logical')
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<1);
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
else
    img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)=...
        (img1(max(1,x-floor(w/2)):min(M,x+floor(w/2)),max(1,y-floor(w/2)):min(N,y+floor(w/2)),:)<128)*255;
    img1(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:)=...
        img(x-floor(w/2)+1:x+floor(w/2)-1,y-floor(w/2)+1:y+floor(w/2)-1,:);
end


function [xs ys W] = smoothing(x,y,L,curve_mode,gau,W);

if L>W
    if curve_mode=='loop' % wrap around the curve by W pixles at both ends
        x1 = [x(L-W+1:L);x;x(1:W)];
        y1 = [y(L-W+1:L);y;y(1:W)];
    else % extend each line curve by W pixels at both ends
        x1 = [ones(W,1)*2*x(1)-x(W+1:-1:2);x;ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
        y1 = [ones(W,1)*2*y(1)-y(W+1:-1:2);y;ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
    end
    
    xx=conv(x1,gau);
    xs=xx(2*W+1:L+2*W);
    yy=conv(y1,gau);
    ys=yy(2*W+1:L+2*W);    
else
    xs = [];
    ys = [];    
end
