function plotStochastic(varargin);
A1 = load (varargin{1});
A2 = load (varargin{2});
A3 = load (varargin{3});
A4 = load (varargin{4});

time11 = A1(:,1);
time21 = A2(:,1);
time31 = A3(:,1);
time41 = A4(:,1);

cons11 = A1(:,2);
cons21 = A2(:,2);
cons31 = A3(:,2);
cons41 = A4(:,2);

time17 = A1(:,3);
time27 = A2(:,3);
time37 = A3(:,3);
time47 = A4(:,3);

cons17 = A1(:,4);
cons27 = A2(:,4);
cons37 = A3(:,4);
cons47 = A4(:,4);

grey=[0.4,0.4,0.4];
plot(time11, cons11, 'r', time21, cons21, 'b', time31, cons31, 'g', time41, cons41, 'b');
