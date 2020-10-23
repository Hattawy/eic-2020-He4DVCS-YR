void Find_CFF( double xB, double t, double &Im, double &Re){
   
 if( 0.1< xB && xB<0.145 ){
          Im  = (99.0 - pow(5.5*t,6)) * exp(-13.47 *t );
          Re  = (-5.6 + 14.7*t - pow(-1.71*t,4)) * exp(-3.9 *t );
          }
    else if( 0.145< xB && xB<0.195 ){
             Im  = (60.77 - pow(5.08 *t, 6)) * exp( -13.2*t );
             Re  = (-8.62 + 15.86*t - pow( -3.59*t, 3)) * exp( -7.27*t);
             }
       else if( 0.195< xB && xB<0.245 ){
                Im  = (39.0 - pow( 5.45*t, 5)) * exp( -12.83*t );
                Re  = (-12.53 + 33.88*t - pow( -0.11*t, 3)) * exp( -7.4*t);
                }
          else if( 0.245< xB && xB<0.295 ){
                   Im  = (26.06 - pow( 4.4*t, 6)) * exp( -12.8*t );
                   Re  = (-18.0 + 35.18*t - pow( -4.56*t, 3)) * exp( -9.5*t);
                   }
             else if( 0.295< xB && xB<0.345 ){
                      Im  = (17.01 - pow( 4.08*t, 6)) * exp( -12.5*t );
                      Re  = (-23.9 + 45.1*t - pow( -5.05*t, 3)) * exp( -10.1*t);
                      }
                else if( 0.345< xB ){
                         Im  = (11.1 - pow( 4.08*t, 6)) * exp( -12.1*t );
                         Re  = (-30.2  -154.7*t - pow( 0.037*t, 3)) * exp( -16.8*t);
                         }
 }
