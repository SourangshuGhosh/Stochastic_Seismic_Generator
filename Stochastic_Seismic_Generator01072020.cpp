#include<iostream>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<stdlib.h>
#include<limits>
#include<vector>
#include<cmath>
#include <sstream> 
#include <fstream>
#include <string.h>
#define two_pi 6.2830
#define pi 3.1415

using namespace std;
string Time,word;
char buf[50];
fstream filenmain2;

void comments_check()
{
	int s,comflag=0;
	filenmain2>>Time;
	stringstream geekwe(Time);
	geekwe>>buf;
	//puts(buf);
	s=strlen(buf);
	
	
			if((buf[0]=='/'&&buf[1]=='*'))
			{
			
				if(buf[s-2]!='*'||buf[s-1]!='/')
				{
			//		puts(buf);
					filenmain2>>Time;
					stringstream geekwqe1(Time);
					geekwqe1>>buf;
			//		puts(buf);
					s=strlen(buf);
				}
				else
				{
			//		puts(buf);
			//		cout<<"Enigma1"<<endl;
					return;
				}
			
				while((buf[s-2]!='*'||buf[s-1]!='/'))
				{
				
					filenmain2>>Time;
					stringstream geekwqe2(Time);
					geekwqe2>>buf;
			//		puts(buf);
					s=strlen(buf);
					//comflag++;
				}
				filenmain2>>Time;
			} 
			//puts(buf);
			//cout<<"Enigma2"<<endl;
}

int main()
{
	// We will first derive the constant Here Reference Distance is set to 1 Km
	// the partition of total shear-wave energy into horizontal components is taken as (1/sqrt(2))
	// F is the effect of the free surface (taken as 2 in almost all applications)
	float fmax=25,fmin=1,cornerfrequency,RadiationPattern,density,Shearwavevelocity,Radius,Constant,Path_effect,Amplification,Dimunution,Site_attenuation,sum;
	float GeometricalSpreading,stressdrop,Momentmagnitude,omega,frequency,TotalSpectralEnergy[30000],Groundmotion_time,Energy,seismicmagnitude,meansquare,rand1[30000];
	float sourcedensity,sourceimpedance,sitedensity,siteimpedance,seismicvelocity,time,groundmotiontime,z0,z1,r1,r2,groundaccln[30000],a,b,c,expandedtime,windowfunction;
	float Seismicsignal[30000],faultarea,faultlength,faultwidth,subfaultwidth,subfaultlength,p1,R1,eta,fTgm,R0_eq7,V_eq7,F_eq7;                                                                                              
	int Numberofsubfaults,NumberLength,NumberWidth,cummulativenumberofsubfaults;
	float QFactor_A,QFactor_B,kappa,Nyquistsamplerate,cornerfrequency1,fouriertransform[30000],deltaf,epsilon,energy[30000];
	float timespentbefore,Kappa,pulsingnoofsubfaults,xdistance,ydistance,zdistance,DipAngle,azimuth,faultstrike,HypocentralDistance,FocalDepth,autocorrelation[30000],time1;
	int pulsecounter;
	char inproblem[30];
	
	cout<<"Please enter the input problemname"<<endl;
		fscanf(stdin,"%s",inproblem);
		
		cout<<"The input filename is"<<endl;
		puts(inproblem);
		
		//filenmain.open(filenameinpproblemname);
		filenmain2.open(inproblem);
		fflush(stdin);
		
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3(Time);
		 getline(geekpo3,word,',');
		 stringstream geek31po(word);
		 geek31po>>RadiationPattern;
			
		 cout<<"\nThe S-wave radiation coefficient is "<<RadiationPattern<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo31(Time);
		 getline(geekpo31,word,',');
		 stringstream geek311po(word);
		 geek311po>>stressdrop;
			
		 cout<<"\nThe difference of stress in the fault plane before and after the rupture occurs is "<<stressdrop<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpod3(Time);
		 getline(geekpod3,word,',');
		 stringstream geek31pod(word);
		 geek31pod>>Shearwavevelocity;
			
		 cout<<"\nThe shear wave velocity(in km/s) is "<<Shearwavevelocity<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpot3(Time);
		 getline(geekpot3,word,',');
		 stringstream geek31epo(word);
		 geek31epo>>density;
			
		 cout<<"\nThe density at the focal point of the earthquake in g/cm^3 is "<<density<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpoq3(Time);
		 getline(geekpoq3,word,',');
		 stringstream geek31poq(word);
		 geek31poq>>Momentmagnitude;
			
		 cout<<"\nThe moment magnitude of the earthquake is "<<Momentmagnitude<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekepo3(Time);
		 getline(geekepo3,word,',');
		 stringstream geeke31po(word);
		 geeke31po>>HypocentralDistance;
			
		 cout<<"\nThe Distance between Epicenter and the site in Km is "<<HypocentralDistance<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3w(Time);
		 getline(geekpo3w,word,',');
		 stringstream geek31pow(word);
		 geek31pow>>FocalDepth;
			
		 cout<<"\nthe Depth of the Focal Point below the earth's crust in Km is "<<FocalDepth<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3qa(Time);
		 getline(geekpo3qa,word,',');
		 stringstream geek31poqa(word);
		 geek31poqa>>DipAngle;
			
		 cout<<"\nThe Dip of the Fault plane in degress is "<<DipAngle<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3sd(Time);
		 getline(geekpo3sd,word,',');
		 stringstream geek31posd(word);
		 geek31posd>>faultstrike;
			
		 cout<<"\nThe Fault Strike of the Fault plane in degrees is "<<faultstrike<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3az(Time);
		 getline(geekpo3az,word,',');
		 stringstream geek31poaz(word);
		 geek31poaz>>azimuth;
			
		 cout<<"\nThe azimuth of the site relative to the origin point(Leftmost corner of the fault rupture square) is "<<azimuth<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b(Time);
		 getline(geekpo3b,word,',');
		 stringstream geek31pob(word);
		 geek31pob>>kappa;
			
		 cout<<"\nThe Dimunution frequency coefficient is "<<kappa<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1(Time);
		 getline(geekpo3b1,word,',');
		 stringstream geek31pob1(word);
		 geek31pob1>>Groundmotion_time;
			
		 cout<<"\nThe Ground Motion Time is "<<Groundmotion_time<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geek7(Time);
		 getline(geek7,word,',');
		 stringstream geek3q(word);
		 geek3q>>QFactor_A;
			
		 cout<<"\nThe Qfactor_A in the form of (QFactor_A)*frequency^(QFactor_B) for the path from source to site is "<<QFactor_A<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geek(Time);
		 getline(geek,word,',');
		 stringstream geek312(word);
		 geek312>>QFactor_B;
			
		 cout<<"\nThe Qfactor_B in the form of (QFactor_A)*frequency^(QFactor_B) for the path from source to site is "<<QFactor_B<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3hj(Time);
		 getline(geekpo3hj,word,',');
		 stringstream geek31pohj(word);
		 geek31pohj>>pulsingnoofsubfaults;
			
		 cout<<"\nThe pulsing percentage rate of the rupture of the fault is "<<pulsingnoofsubfaults<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3ba1(Time);
		 getline(geekpo3ba1,word,',');
		 stringstream geekr31pob1(word);
		 geekr31pob1>>fmin;
			
		 cout<<"\nThe minimum frequecny of the earthquake is "<<fmin<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpoa3b14(Time);
		 getline(geekpoa3b14,word,',');
		 stringstream geeka31pob1z(word);
		 geeka31pob1z>>fmax;
			
		 cout<<"\nThe maximum frequecny of the earthquake is "<<fmax<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpoa3b1(Time);
		 getline(geekpoa3b1,word,',');
		 stringstream geek31apob1(word);
		 geek31apob1>>deltaf;
			
		 cout<<"\nThe value of the frequency increment is "<<deltaf<<endl;
		 
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1z(Time);
		 getline(geekpo3b1z,word,',');
		 stringstream geek31pob1z(word);
		 geek31pob1z>>epsilon;
			
		 cout<<"\nThe Value of Window Parameter epsilon is "<<epsilon<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1za(Time);
		 getline(geekpo3b1za,word,',');
		 stringstream geek31pob1za(word);
		 geek31pob1za>>eta;
			
		 cout<<"\nThe Value of Window Parameter eta is "<<eta<<endl;
		 
		  comments_check();
		 comments_check();
		 
		 stringstream geekpo3sb1(Time);
		 getline(geekpo3sb1,word,',');
		 stringstream geek31spob1(word);
		 geek31spob1>>fTgm;
			
		 cout<<"\nThe Value of Window Parameter fTgm is "<<fTgm<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1ax(Time);
		 getline(geekpo3b1ax,word,',');
		 stringstream geek31pob1zx(word);
		 geek31pob1zx>>R0_eq7;
			
		 cout<<"\nThe Reference distance of the source is "<<R0_eq7<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1g(Time);
		 getline(geekpo3b1g,word,',');
		 stringstream geek31pob1g(word);
		 geek31pob1g>>R1;
			
		 cout<<"\nThe Value of Geometrical Spreading parameter 1 is "<<R1<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1q(Time);
		 getline(geekpo3b1q,word,',');
		 stringstream geek31pob1ag(word);
		 geek31pob1ag>>p1;
			
		 cout<<"\nThe Value of Geometrical Spreading parameter 2 is "<<p1<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1qj(Time);
		 getline(geekpo3b1qj,word,',');
		 stringstream geek31pob1agsq(word);
		 geek31pob1agsq>>V_eq7;
			
		 cout<<"\nThe partition of total shear-wave energy into horizontal components is "<<V_eq7<<endl;
		 
		 comments_check();
		 comments_check();
		 
		 stringstream geekpo3b1q12(Time);
		 getline(geekpo3b1q12,word,',');
		 stringstream geek31pob1agwe(word);
		 geek31pob1agwe>>F_eq7;
			
		 cout<<"\nThe effect of the free surface is "<<F_eq7<<endl;
		 //cout<<"cos(90)"<<cos(90)<<endl;
		 //all input parameters are read 
		 
		 	
	Constant=RadiationPattern*V_eq7*F_eq7/(4*pi*density*pow(Shearwavevelocity,3)*R0_eq7);
	cout<<"Constant:"<<Constant<<endl; // Eq (7) in Boore 2003
	
	seismicmagnitude=pow(10,(1.5*(Momentmagnitude+10.7))); // M0 in Eq 2 of Boore 2003
	cout<<"Seismic_Magnitude:"<<seismicmagnitude<<endl;
	//calculating the subfault dimension by Wells and Coppersmith Formula(1994)
	faultarea=exp(-3.49+0.91*Momentmagnitude);
	faultlength=exp(-1.88+0.5*Momentmagnitude);	
	faultwidth=faultarea/faultlength;
	cout<<"Faultlength:"<<faultlength<<endl;
	cout<<"Faultwidth"<<faultwidth<<endl;
	
	//calculating the subfault dimension by Beresnev and Atkinson(2002)
	subfaultwidth=-2+0.4*Momentmagnitude;
	subfaultlength=-2+0.4*Momentmagnitude;
	cout<<"SubFaultlength:"<<subfaultlength<<endl;
	cout<<"SubFaultwidth"<<subfaultwidth<<endl;
		
	Numberofsubfaults=faultarea/(subfaultwidth*subfaultlength);// Boore (2009)
	NumberLength=faultlength/subfaultlength;
	NumberWidth=faultwidth/subfaultwidth;
	pulsingnoofsubfaults=(pulsingnoofsubfaults/100)*Numberofsubfaults;
	pulsingnoofsubfaults=(int)pulsingnoofsubfaults;

	cout<<"Numberofsubfaults:"<<Numberofsubfaults<<endl;
	cout<<"Number of subfault units along the Length is "<<NumberLength<<endl;
	cout<<"Number of subfault units along the Width is "<<NumberWidth<<endl;
	cout<<"Pulsing No. of Subfaults is "<<pulsingnoofsubfaults<<endl;
	
	cummulativenumberofsubfaults=pulsingnoofsubfaults;
		
	timespentbefore=0;
	pulsecounter=0;
		
	int counter=0;
	for(frequency=fmin;frequency<fmax;frequency=frequency+deltaf) //in hertz 
	{
			TotalSpectralEnergy[counter]=0;	
			counter++;
	}
	
	for(int i=0;i<NumberLength;i++)
	{
		for(int j=0;j<NumberWidth;j++)
		{
			cout<<"Finding the contribution of earthquake from subfault at position ["<<i<<"]["<<j<<"]"<<endl;
			cout<<"Cummulative number of subfaults:"<<cummulativenumberofsubfaults<<endl;
			cout<<"Pulsecounter:"<<pulsecounter<<endl;
			cout<<"Pulsingnoofsubfaults"<<pulsingnoofsubfaults<<endl;
			if(pulsecounter==pulsingnoofsubfaults)  // if pulsecounter becomes equal to the pulsingnoofsubfaults( it means the maximum number of subfaults that can be active together)
			{
				pulsecounter=0;
				cummulativenumberofsubfaults=cummulativenumberofsubfaults+pulsingnoofsubfaults;
			}
			
			//calculating the distance of the subfault center to the site
			xdistance=HypocentralDistance*cos(((azimuth-faultstrike)/180)*pi)-(2*j-1)*(subfaultlength/2);
			ydistance=HypocentralDistance*sin(((azimuth-faultstrike)/180)*pi)-(2*i-1)*(subfaultwidth/2)*sin(((DipAngle)/180)*pi);
			zdistance=FocalDepth+(2*i-1)*(subfaultwidth/2)*cos(((DipAngle)/180)*pi);
			
			Radius=sqrt(pow(xdistance,2)+pow(ydistance,2)+pow(zdistance,2));
			cout<<"Total Distance of this subfault from the site is "<<Radius<<endl;
			
			cornerfrequency=4.9E06*pow((Numberofsubfaults/cummulativenumberofsubfaults),0.33333333)*Shearwavevelocity*pow(((Numberofsubfaults*stressdrop)/(seismicmagnitude)),0.33333333)*10;//*pow(((Numberofsubfaults*Numberofsubfaults)/(cummulativenumberofsubfaults)),(1/3));
			    // Eq 4 in Boore 2003
			cout<<"The corner frequency of the subfault is "<<cornerfrequency<<endl;
			
			if(Radius<R1)
			{
				GeometricalSpreading=R0_eq7/Radius; // Eq 8 of Boore 2003
			}
			else
			{
				GeometricalSpreading=R0_eq7/R1*pow((R1/Radius),p1);
			}
			cout<<"GeometricalSpreading:"<<GeometricalSpreading<<endl;
			cout<<"\n\n"<<endl;
			// We first calculate the spectral energy content of the earthquake wave near this subfault whose index position is i and j
			counter=0;
			
			for(frequency=fmin;frequency<fmax;frequency=frequency+deltaf)
			{
				Energy=(Constant*(seismicmagnitude/Numberofsubfaults)*1)/(1+pow((frequency/cornerfrequency),2))*pow(10,-20); // Constant defined above
				//cout<<"Energy at frquency:"<<frequency<<" is "<<Energy<<endl;
			
				Path_effect=GeometricalSpreading*exp(-(pi*frequency*Radius)/(Shearwavevelocity*QFactor_A*pow(frequency,QFactor_B))); // Eq 8 in Boore 2003
			//	cout<<"Path_Effect:"<<Path_effect<<endl;
				
				if(frequency<=0.5) // Eq 10
				{
					Amplification=1;
				}
				else if(frequency<=1)
				{
					Amplification=(1.5-1)*2*(frequency-0.5)+1;
				}
				else if(frequency<=2)
				{
					Amplification=(2-1.5)*1*(frequency-1)+1.5;
				}
				else if(frequency<=5)
				{
					Amplification=(2.5-2)*(1/3)*(frequency-2)+2;
				}
				else if(frequency<10)
				{
					Amplification=(3-2.5)*(1/5)*(frequency-5)+2.5;
				}
				else
				{
					Amplification=(3-2.5)*(1/5)*(frequency-5)+2.5;
				}
				
				Dimunution=exp(-(two_pi/2)*kappa*frequency); // Eq 20
				
				Site_attenuation=Amplification*Dimunution;
				//cout<<"Site_Attenuation"<<Site_attenuation<<endl;
				
				omega=two_pi*frequency;
				TotalSpectralEnergy[counter]=TotalSpectralEnergy[counter]+pow((omega*omega*Energy*Path_effect*Site_attenuation*((1+(frequency/cornerfrequency1)*(frequency/cornerfrequency1))/(1+(frequency/cornerfrequency1)*(frequency/cornerfrequency1)))*deltaf*2*pi*sqrt(Numberofsubfaults)),2);      
				//cout<<"Total Spectral Energy:"<<TotalSpectralEnergy[counter]<<endl;
				counter++;
			}
			pulsecounter++;
		}
	}
						
			
			//If a function x(t) contains no frequencies higher than B hertz, it is completely determined by giving its ordinates at a series of points spaced 1/2B seconds apart.
			Nyquistsamplerate=(1/(2*fmax));
			counter=0;
			// We fnow generate a white noise with 0 mean and standard devation of 1
			
			FILE *fptro1=fopen("Earthquake1.csv","w");
			for(time=0;time<Groundmotion_time;time=time+Nyquistsamplerate)
			{
				r1= ((double) rand() / (RAND_MAX));
				r2= ((double) rand() / (RAND_MAX));
									
				z0 = sqrt(-2.0 * log(r1)) * cos(two_pi * r2);       //Using Box-Mueller Theorem to generate two independent normal variables having 0 mean
				z1 = sqrt(-2.0 * log(r1)) * sin(two_pi * r2);
				groundaccln[counter]=z0;
				
				fprintf(fptro1,"%6.9f,%6.9f\n",time,groundaccln[counter]);
				counter++;	
			}
			fclose(fptro1);
			//Now we convert the white noise to a windowed noise by multiplying a window function
			
			meansquare=0;
			int totalcount=0;
			b=-((epsilon*log(eta))/(1+epsilon*(log(epsilon)-1)));
			c=b/epsilon;
			a=pow((exp(1)/eta),b);
			expandedtime=fTgm*Groundmotion_time; // Eq28 of Boore 2003
			
			cout<<"a:"<<a<<endl;
			cout<<"b:"<<b<<endl;
			cout<<"c:"<<c<<endl;
			cout<<"expandedtime"<<expandedtime<<endl;
			
			FILE *fptro2=fopen("Earthquake2.csv","w");
			for(time=0;time<Groundmotion_time;time=time+Nyquistsamplerate)
			{
				windowfunction=a*pow((time/expandedtime),b)*exp(-c*(time/expandedtime));
				groundaccln[totalcount]=windowfunction*groundaccln[totalcount];
				
				if(isinf(groundaccln[totalcount]))
				{
					groundaccln[totalcount]=0;
				}
				fprintf(fptro2,"%6.9f,%6.9f\n",time,groundaccln[totalcount]);
				cout<<groundaccln[totalcount]<<endl;
				meansquare=meansquare+groundaccln[totalcount]*groundaccln[totalcount];
				totalcount++;
			}
			fclose(fptro2);	
			cout<<totalcount<<endl;
			
			meansquare=meansquare/totalcount;
			meansquare=sqrt(meansquare);
			cout<<"The meansqare is"<<meansquare<<endl;
			counter=0;
			
			
			//converting the signal back to frequency domain
			
			//first finding the autocorrelation function
			
			counter=0;
			
			for(time=0;time<Groundmotion_time;time=time+Nyquistsamplerate) // for any lag tau, please divide the sum by N, not N-tau
			{
				totalcount=0;
				autocorrelation[counter]=0;
				for(time1=0;time1<(Groundmotion_time-time);time1=time1+Nyquistsamplerate)
				{
					autocorrelation[counter]=autocorrelation[counter]+groundaccln[totalcount]*groundaccln[totalcount+counter];
					totalcount++;
				}
				autocorrelation[counter]=autocorrelation[counter]/(totalcount+(time/Nyquistsamplerate));
				counter++;
			}
			
			char jk;
			//Now finding out the spectrum of the ACF of the windowed function
			counter=0;
			float deltaomega=two_pi*deltaf;
			cout<<deltaomega<<endl;
			
			FILE *fptro3=fopen("Earthquake3.csv","w");
			FILE *fptro4=fopen("Earthquake4.csv","w");
			FILE *fptro5=fopen("Earthquake5.csv","w");
			
			for(frequency=fmin;frequency<fmax;frequency=frequency+deltaf)
			{
				fouriertransform[counter]=0;
				totalcount=0;
				for(time=0;time<Groundmotion_time;time=time+Nyquistsamplerate)
				{
					fouriertransform[counter]=fouriertransform[counter]+autocorrelation[totalcount]*cos(two_pi*frequency*time);
					totalcount++;
				}
				fouriertransform[counter]=(1/(two_pi))*fouriertransform[counter]*Nyquistsamplerate;
				
				if(fouriertransform[counter]<0)
				fouriertransform[counter]=fouriertransform[counter-1]; //why?
				
				TotalSpectralEnergy[counter]=sqrt(TotalSpectralEnergy[counter]/Numberofsubfaults);
				
				energy[counter]=(TotalSpectralEnergy[counter]*fouriertransform[counter]/meansquare);
				
				fprintf(fptro3,"%6.9f,%6.9f\n",frequency,fouriertransform[counter]);
				fprintf(fptro4,"%6.9f,%6.9f\n",frequency,energy[counter]);
				fprintf(fptro5,"%6.9f,%6.9f\n",frequency,TotalSpectralEnergy[counter]);
				
				cout<<fouriertransform[counter]<<","<<TotalSpectralEnergy[counter]<<","<<energy[counter]<<endl;
				counter++;
			}
			fclose(fptro3);
			fclose(fptro4);
			fclose(fptro5);
			
			cin>>jk;
			//multiply with ground motion spectrum
			counter=0;
			for(frequency=fmin;frequency<fmax;frequency=frequency+deltaf)
			{
				r1= ((double) rand() / (RAND_MAX));
				rand1[counter]=r1;	
				counter++;
			}
			int counter2=0;
			
			
			FILE *fptro=fopen("Earthquake.csv","w");
			for(time=0;time<Groundmotion_time;time=time+Nyquistsamplerate)
			{
				sum=0;
				counter=0;
				for(frequency=fmin;frequency<fmax;frequency=frequency+deltaf)
				{
					//cout<<"Energy"<<energy[counter]<<endl;
					sum=sum+sqrt(2*energy[counter])*cos(two_pi*frequency*time+two_pi*rand1[counter]);	
					counter++;
				}
				Seismicsignal[counter2]=sum;
						
				fprintf(fptro,"%6.9f,%6.9f\n",time,Seismicsignal[counter2]);
				cout<<Seismicsignal[counter2]<<endl;
				counter2++;
			}
			fclose(fptro);
	
	return 0;
}
