#ifndef SKELANA_PSCTBL_HPP
#define SKELANA_PSCTBL_HPP

namespace skelana
{
    inline const int NVMAX = 511;
    extern "C" struct
    {
      int npa;                        
      int nst;                       
      int nsh;                     
      int nlu;                   
      
      int ishst[NVMAX];
      int istsh[NVMAX];

      int istpa[NVMAX];
      int ipast[NVMAX];
      
      int ishlu[NVMAX];
      int ilush[NVMAX];

      int ilust[NVMAX];
      int istlu[NVMAX];

      int ipapv[NVMAX][2];
      int istvx[NVMAX][2];
      
    } psctbl_;

  inline int& NPA () { return psctbl_.npa; }
  inline int& NST () { return psctbl_.nst; }
  inline int& NSH () { return psctbl_.nsh; }
  inline int& NLU () { return psctbl_.nlu; }
  
  inline int& ISHST(int j) { return psctbl_.ishst[j-1]; } 
  inline int& ISTSH(int j) { return psctbl_.istsh[j-1]; }
  
  inline int& ISTPA(int j) { return psctbl_.istpa[j-1]; }
  inline int& IPAST(int j) { return psctbl_.ipast[j-1]; }

  inline int& ISHLU(int j) { return psctbl_.ishlu[j-1]; }
  inline int& ILUSH(int j) { return psctbl_.ilush[j-1]; }

  inline int& ILUST(int j) { return psctbl_.ilust[j-1]; }
  inline int& ISTLU(int j) { return psctbl_.istlu[j-1]; }

  inline int& IPAPV(int i,int j) { return psctbl_.ipapv[j-1][i-1]; }
  inline int& ISTVX(int i,int j) { return psctbl_.istvx[j-1][i-1]; }

#endif // SKELANA_PSCTBL_HPP
