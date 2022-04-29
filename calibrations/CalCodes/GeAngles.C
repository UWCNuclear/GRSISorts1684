// Printing angles from GRSISORT/librairies/TGRSIAnalysis/TTigress/TTigress.cxx and /TGriffin/TGriffin.cxx
// Note that TIGRESS lampshades are rotated 45 degrees from GRIFFIN!!!!

{
		printf("GRIFFIN\n");
  		for(int k=1;k<=16;k++) {

// TVector3 TGriffin::GetPosition(int DetNbr, int CryNbr, double dist)
	          Double_t dtheta = TGriffin::GetPosition(k).Theta()*TMath::RadToDeg();
		  Double_t dphi = TGriffin::GetPosition(k).Phi()*TMath::RadToDeg();
 
  		  printf("Detector %i :   Theta =   %8.2f      Phi =   %8.2f   \n\n\n", k, dtheta, dphi);
 
 		}

		printf("GRIFFIN 11 cm\n");
  		for(int i=1;i<=16;i++) {
   		  for(int j=0;j<4;j++) {

     		  int det1 = i*4+j;
// TVector3 TTigress::GetPosition(int DetNbr, int CryNbr, int SegNbr, double dist, bool smear)
	          Double_t gtheta = TGriffin::GetPosition(i,j,110.).Theta()*TMath::RadToDeg();
		  Double_t gphi = TGriffin::GetPosition(i,j,110.).Phi()*TMath::RadToDeg();
 
  		  printf("Detector %i, Crystal %i :\t\t  Theta =   %8.2f\t\t      Phi =   %8.2f   \n\n\n", i, det1-3, gtheta, gphi);
    
    		  }
 		}

		printf("TIGRESS 14.5 cm\n");
  		for(int i=1;i<=16;i++) {
   		  for(int j=0;j<4;j++) {

     		  int det1 = i*4+j;
// TVector3 TTigress::GetPosition(int DetNbr, int CryNbr, int SegNbr, double dist, bool smear)
	          Double_t ctheta = TTigress::GetPosition(i,j,0,145.).Theta()*TMath::RadToDeg();
		  Double_t cphi = TTigress::GetPosition(i,j,0,145.).Phi()*TMath::RadToDeg();
 
  		  printf("Detector %i, Crystal %i :\t\t   Theta =   %8.2f\t\t      Phi =   %8.2f   \n", i, det1-3, ctheta, cphi);
    
    		  }
 		}
}
