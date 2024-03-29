Line.con = [ ... 
  2  5  100  69  60  0  0  0.05695  0.17388  0.034  0  0  0  0  0  1;
  6  12  100  13.8  60  0  0  0.12291  0.25581  0  0  0  0  0  0  1;
  12  13  100  13.8  60  0  0  0.22092  0.19988  0  0  0  0  0  0  1;
  6  13  100  13.8  60  0  0  0.06615  0.13027  0  0  0  0  0  0  1;
  6  11  100  13.8  60  0  0  0.09498  0.1989  0  0  0  0  0  0  1;
  11  10  100  13.8  60  0  0  0.08205  0.19207  0  0  0  0  0  0  1;
  9  10  100  13.8  60  0  0  0.03181  0.0845  0  0  0  0  0  0  1;
  9  14  100  13.8  60  0  0  0.12711  0.27038  0  0  0  0  0  0  1;
  14  13  100  13.8  60  0  0  0.17093  0.34802  0  0  0  0  0  0  1;
  7  9  100  13.8  60  0  0  0  0.11001  0  0  0  0  0  0  1;
  1  2  100  69  60  0  0  0.01938  0.05917  0.0528  0  0  0  0  0  1;
  3  2  100  69  60  0  0  0.04699  0.19797  0.0438  0  0  0  0  0  1;
  3  4  100  69  60  0  0  0.06701  0.17103  0.0346  0  0  0  0  0  1;
  1  5  100  69  60  0  0  0.05403  0.22304  0.0492  0  0  0  0  0  1;
  5  4  100  69  60  0  0  0.01335  0.04211  0.0128  0  0  0  0  0  1;
  2  4  100  69  60  0  0  0.05811  0.17632  0.0374  0  0  0  0  0  1;
  5  6  100  69  60  0  5  0  0.25202  0  0.932  0  0  0  0  1;
  4  9  100  69  60  0  5  0  0.55618  0  0.969  0  0  0  0  1;
  4  7  100  69  60  0  5  0  0.20912  0  0.978  0  0  0  0  1;
  8  7  100  18  60  0  1.304348  0  0.17615  0  0  0  0  0  0  1;
 ];

Bus.con = [ ... 
  1  69  1  0  4  1;
  2  69  1  0  4  1;
  3  69  1  0  4  1;
  4  69  1  0  4  1;
  5  69  1  0  4  1;
  6  13.8  1  0  2  1;
  7  13.8  1  0  2  1;
  8  18  1  0  3  1;
  9  13.8  1  0  2  1;
  10  13.8  1  0  2  1;
  11  13.8  1  0  2  1;
  12  13.8  1  0  2  1;
  13  13.8  1  0  2  1;
  14  13.8  1  0  2  1;
 ];

SW.con = [ ... 
  1  100  69  1.06  0  9.9  -9.9  1.2  0.8  2.324  1  1  1;
 ];

PV.con = [ ... 
  2  100  69  0.4  1.045  0.5  -0.4  1.2  0.8  1  1;
  6  100  13.8  0  1.07  0.24  -0.06  1.2  0.8  1  1;
  3  100  69  0  1.01  0.4  0  1.2  0.8  1  1;
  8  100  18  0  1.09  0.24  -0.06  1.2  0.8  1  1;
 ];

PQ.con = [ ... 
  11  100  13.8  0.049  0.0252  1.2  0.8  0  1;
  13  100  13.8  0.189  0.0812  1.2  0.8  0  1;
  3  100  69  1.3188  0.266  1.2  0.8  0  1;
  5  100  69  0.1064  0.0224  1.2  0.8  0  1;
  2  100  69  0.3038  0.1778  1.2  0.8  0  1;
  6  100  13.8  0.1568  0.105  1.2  0.8  0  1;
  4  100  69  0.6692  0.056  1.2  0.8  0  1;
  14  100  13.8  0.2086  0.07  1.2  0.8  0  1;
  12  100  13.8  0.0854  0.0224  1.2  0.8  0  1;
  10  100  13.8  0.126  0.0812  1.2  0.8  0  1;
  9  100  13.8  0.413  0.2324  1.2  0.8  0  1;
 ];

Syn.con = [ ... 
  1  615  69  60  5.2  0.2396  0  0.8979  0.6  0.23  7.4  0.03  0.646  0.646  0.4  0  0.033  10.296  2  0  0  1  1  0  0  0;
  3  60  69  60  6  0  0.0031  1.05  0.185  0.13  6.1  0.04  0.98  0.36  0.13  0.3  0.099  13.08  2  0  0  1  1  0  0  0;
  2  60  69  60  6  0  0.0031  1.05  0.185  0.13  6.1  0.04  0.98  0.36  0.13  0.3  0.099  13.08  2  0  0  1  1  0  0  0;
  8  25  18  60  6  0.134  0.0014  1.25  0.232  0.12  4.75  0.06  1.22  0.715  0.12  1.5  0.21  10.12  2  0  0  1  1  0  0  0;
  6  25  13.8  60  6  0.134  0.0014  1.25  0.232  0.12  4.75  0.06  1.22  0.715  0.12  1.5  0.21  10.12  2  0  0  1  1  0  0  0;
 ];

Exc.con = [ ... 
  1  2  7.32  0  200  0.02  0.002  1  1  0.2  0.001  0.0006  0.9;
  3  2  4.38  0  20  0.02  0.001  1  1  1.98  0.001  0.0006  0.9;
  2  2  4.38  0  20  0.02  0.001  1  1  1.98  0.001  0.0006  0.9;
  4  2  6.81  1.395  20  0.02  0.001  1  1  0.7  0.001  0.0006  0.9;
  5  2  6.81  1.395  20  0.02  0.001  1  1  0.7  0.001  0.0006  0.9;
 ];

Bus.names = {... 
  'Bus 01'; 'Bus 02'; 'Bus 03'; 'Bus 04'; 'Bus 05'; 
  'Bus 06'; 'Bus 07'; 'Bus 08'; 'Bus 09'; 'Bus 10'; 
  'Bus 11'; 'Bus 12'; 'Bus 13'; 'Bus 14'};

