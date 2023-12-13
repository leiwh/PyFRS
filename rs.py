
#================================================
#==============Reverse shock=====================

#-----------Scientific Library-------
import math
import numpy as np
import scipy.optimize

#-----------My Library---------
#Credit: Weihua Lei
import cgs
import astro
import grb


#def RS_flux(self):
def RS_flux(time_obs,nu_obs,**Z):
    log_obs_time = np.log10(time_obs)

#-------Input parameters: independent of model
#        zi= z
#        D28=astro.LUD(zi)/1.e28

#-------Unpack parameters
    for key,value in Z.items():
        if key == 'z':              # Redshift
            z = value
        elif key == 'k':            # Medium density profile index
            k = value
        #elif key == 'SSC':          # SSC
        #    if value == 1:
        #        SSC = 'Yes'
        #    else:
        #        SSC = 'No'
        #elif key == 'XIC':          # EIC
        #    if value == 1:
        #        XIC = 'Yes'
        #    else:
        #        XIC = 'No'
        elif key == 'E0':           # Isotropic-equivalent energy in 10^52 erg
            E52 = value/1.0e52
        elif key == 'Gamma0':       # Initial bulk Lorentz factor
            Gm0 = value
        elif key == 'theta_j':      # Jet opening angle
            theta_j = value
        elif key == 'theta_obs':    # Viewing angle in deg
            theta_obs = value
        elif key == 'n18':          # CNM density at 10^18cm
            n18 = value
        elif key == 'p':            # Electron energy distribution index
            p = value
        elif key == 'epsilon_e':    # Fraction of the shock energy into electrons
            epsilon_e = value
        elif key == 'epsilon_B':    # Fraction of the shock energy in magnetic field
            epsilon_B = value




#-------Input for jet
    zi=z
    D28=astro.LUD(zi)/1.e28

    n1=n18
    pp=p


    mytime=1.*cgs.day
    


    # Band={'X-ray': 1e18,'Optical': 1e15,'Radio': 1e9, 'GeV': 2.4e23}
#        Time={'input':tday, 'adjust':10.0**lgtday}
#    Dynmodel={'analytical': 1,'differential': 2}
#    with_SSC={'Yes':1, 'No':0}
#    with_jetbreak={'Yes':1, 'No':0}
    Tunits = {'Second': 1.0, 'Day': cgs.day}
    Nuunits= {'Hz': 1.0, 'keV': cgs.keV2Hz}

#    dynamic_model='differential'
    # Band_mode='X-ray'
#        v=Band[Band_mode]
    vlist=list(dict.fromkeys(list(nu_obs)))
    Funit=1.
    time_unit='Second'
    nu_unit='Hz'
#    jet_break='Yes'
#    inj_model='No'
#    Smooth='No'



#-------Initialize the time, gamma, Fv arries in Log space
    inum=500

    t=np.logspace(-4., 3.0,inum)
    to=np.logspace(-4., 3.0,inum)  # t_observe

    Rt=np.logspace(14.,21.,inum)
    Gmt=np.logspace(2.,-5.,inum)

    para=np.zeros((inum,19))

    lgt=np.logspace(-4., 3.0,inum)
    lgFv=np.logspace(1.,6.,inum)*np.ones((len(vlist),inum))
#    lgFv_SSC=np.logspace(1.,6.,inum)*np.ones((len(vlist),inum))

#    lgFv=np.logspace(1.,6.,inum)
#        lgFv_SSC=np.logspace(1.,6.,inum)

#    lgNu=np.logspace(7.,30.,inumv)
#    lgFnu=np.logspace(1.,6.,inumv)
#        lgFnu_SSC=np.logspace(1.,6.,inumv)

#-------Initialize critical parameters
#    lgNua=np.logspace(-4.,3.0,inump)
#    lgNum=np.logspace(-4.,3.0,inump)
#    lgNuc=np.logspace(-4.,3.0,inump)
#    lgFluxa=np.logspace(-4.,3.0,inump)
#    lgFluxm=np.logspace(-4.,3.0,inump)
#    lgFluxc=np.logspace(-4.,3.0,inump)

#    lgtdec=np.logspace(-4.,3.0,inump)
#    lgFluxdec=np.logspace(-4.,3.0,inump)
#         lgtSed=np.logspace(-4.,3.0,inump)
#         lgFluxSed=np.logspace(-4.,3.0,inump)


#-------start calculation
#    Gm0= Gamma0
    Gmt[1]=Gm0

    R0=Rt[1]
    Eiso=1.e52 *E52
#    t90=T90
    t[0]=0.
    t[1]=R0/(2.*Gm0**2.0*cgs.c*cgs.day)
    to[0]=0.
    to[1]=t[1]

    #if ( lgtday<math.log(t[1],10)):
    #     lgtday=math.log(t[1],10)

    Gm3x=Gm0

    An0=n18 *1.e18**k
    lsd=((3.- k)*Eiso/(4.*cgs.pi*An0*cgs.mp*cgs.c**2.))**(1./(3. - k))  #Sedov length
    tx=lsd/(2.*cgs.c* Gm0**(2.+2/(3.- k) ) )    #crossing time, also decelation time
    txo=tx*(1+zi)
    txo0=txo
    Rx=(2.*Gm0**2.0 *cgs.c)*tx
    Ne0=Eiso/(Gm0*cgs.mp*cgs.c**2.)
        
    t[1]=Rt[1]/(2.*Gm0**2.0*cgs.c) *(1.+zi)

#        txo= t_dec_n*Tunits[ time_unit]
    t_dec=txo/Tunits[time_unit]
    tx=txo/(1+zi)
    Rx=(2.*Gm0**2.0 *cgs.c)*tx
        

        
    if (k==0.):
        g=2.
    else:
        g=1.

    n1x=n18*(Rx/1.e18)**(-k)
    n3x=7.*n1x*(lsd/Rx)**(3-k) 

    e3x=4.*Gm0**2. *n1x* cgs.mp* cgs.c**2.
    Bcx=(8.*cgs.pi*e3x*  epsilon_B)**0.5
    fn41x=(lsd/Rx)**(3.- k)
    Gm341x=4.*Gm0**2 *fn41x**(-1)/7.
    Ne3x=Ne0 *Gm0**((3- k)/3.)* (Rx/lsd)**((3- k)/2.)

#-----------Minimum Lorentz factor-------
    gm_mx=Gm341x*  epsilon_e* ( p-2.)/( p-1.) *cgs.mp/cgs.me
#-----------Cooling Lorentez factor
    gm_cx=6.*cgs.pi*cgs.me*cgs.c/(cgs.sigmaT* Bcx**2 * Gm0* txo/(1.+zi))
        

    Pvmx=grb.Pvmax(Gm3x,Bcx)

    #print("### RS:: tx0,gm_cx0,gm_mx0,Rx,n1x0,n3x,e3x,Bcx,fn41x,Gm341x:", txo0,gm_cx,gm_mx,Rx,n1x,n3x,e3x,Bcx,fn41x,Gm341x)

    Gm3=Gm0
    for i in range(1,inum):
        n1= n18*(Rt[i]/1.e18)**(- k)

        fview=1.

        Ri=Rt[i]
        t[i]=Rt[i]/(2.*Gm3**2.0 *cgs.c) *(1.+zi)


#            e3=4.*cgs.Gm0**2. *n1*cgs.mp* cgs.c**2.
        if (t[i]<txo):
            Gm3=Gm0
#                t[i]=Rt[i]/(2.*Gm0**2.0 *cgs.c) *(1.+zi)

            Ne3=Ne0 *Gm0* (Rt[i]/lsd)**((3- k)/2.)
#                Ne3=Ne0* (t[i]/txo)**((3- k)/2.)
#-----------Internal energy---------
            e3=4.*Gm0**2. *n1*cgs.mp* cgs.c**2.
            n3=7.*n1*(lsd/Rt[i])**(3- k) 

#-----------B field at comoving frame-------
#            Bc=grb.Bco2( epsilon_B,e3)
            Bc=(8.*cgs.pi*e3*  epsilon_B)**0.5
            
            fn41=(lsd/Rt[i])**(3.- k)
            Gm341=4.*Gm0**2 *fn41**(-1)/7.

#-----------Minimum Lorentz factor-------
#            gm_m=grb.gamma_m2( epsilon_e,GM21,pp,gm_Max)
            gm_m=Gm341*  epsilon_e* ( p-2.)/( p-1.) *cgs.mp/cgs.me
#                gm_m= epsilon_e* ( p-2.)/( p-1.) *(e3/n3)/(cgs.me*cgs.c**2.)
#-----------Cooling Lorentez factor
#            gm_c=grb.gamma_c(Bc,GM21,ti*fview,zi)
#            gm_c=grb.gamma_c(Bc,Gm0,t[i],zi)
            gm_c=6.*cgs.pi*cgs.me*cgs.c/(cgs.sigmaT* Bc**2* Gm0* t[i]/(1.+zi))

#                Pvm=grb.Pvmax(Gm3,Bc)


#                dL=Ne3/(4.*cgs.pi*Rt[i]**2. *n3)

#                n3x=n3
#                e3x=e3
#                Bcx=Bc
#                gm_mx=gm_m
#                gm_cx=gm_c
#                Rx=Rt[i]
                
        else:
#                t[i]= txo*  (Rt[i]/Rx)**(1+2.*g)
#                t[i]=Rt[i]/(2.*Gm3**2.0*cgs.c) *(1.+zi)
            Gm3=Gm3x* (t[i]/txo)**(-g/(1+2.*g))

            Ne3= Ne0
            n3= n3x *(t[i]/txo)**(-6.*(3+g)/(7 *(1+2.*g)) )
            e3= e3x* (t[i]/txo)**(-8.*(3+g)/(7.*(1+2.*g)) )

#-----------B field at comoving frame-------
#                Bc=(8.*cgs.pi *e3 * epsilon_B)**0.5
            Bc=Bcx* (t[i]/txo)**(-4.*(3+g)/(7.*(1+2.*g)) )
            
            fn41=(lsd/Rt[i])**(3.- k)
            Gm341=4.*Gm0**2 *fn41**(-1)/7.

#-----------Minimum Lorentz factor-------
            gm_m=gm_mx *(t[i]/txo)**(-2.*(3.+g)/(7.*(1+2.*g)) )
#-----------Cooling Lorentez factor
            gm_c=gm_cx *(t[i]/txo)**(-2.*(3.+g)/(7.*(1+2.*g)) )

#                Pvm=Pvmx* (t[i]/txo)**(-4.*(3+g)/(7.*(1+2.*g)) -g/(1+2.*g) )

#                t[i]=Rt[i]/(2.*Gm3**2.0*cgs.c) *(1.+zi)
#                dL=Ne3/(4.*cgs.pi*Rx**2. *n3x)

 

            

        ti=t[i]/fview  #off-viewer's time
        to[i]=ti
        lgt[i]=math.log(ti/Tunits[ time_unit],10)

        ei=e3

        beta=0.1
#            beta=math.sqrt(2.*(Gm3-1))
#            beta=math.sqrt(1.-1/Gm3**2.0)
            
#-----------Max Lorentz factor-------
        gm_Max = grb.gamma_Max(Bc)
             
        
#-----------v_Max, v_m, v_c, v_a----------
        vMax = grb.nu_gme(Gm3,gm_Max, Bc,zi)
        vm=grb.nu_gme(Gm3,gm_m, Bc,zi)
        vc=grb.nu_gme(Gm3,gm_c, Bc,zi)

        dL=Ne3/(4.*cgs.pi*Rt[i]**2. *n3)
        va=grb.nu_aB(Bc,pp,n3/(4.*Gm3),vm,vc,Gm3,dL*Gm3,zi)

        gm_a=(va*(2.*cgs.pi *cgs.me *cgs.c)*(1.+zi)/(Gm3 *cgs.qe *Bc))**0.5
#-----------Y-parameter for SSC without K-N correction
        UB = Bc**2. /(8.*cgs.pi)
        Uph=1.


#-----------t_min, t_col and t_colmin, t_jet
        # if (i>1):
        #     if ((v - vm) * (v - vm0) <= 0.0 ):
        #         #print("### RS:: time for vm cross:", Rt[i-1], ti/Tunits[ time_unit])
        #         t_min=ti/Tunits[time_unit]
        #     if ((v - vc) * (v - vc0) <= 0.0 ):
        #         #print("### RS:: time for vc cross:", Rt[i-1], ti/Tunits[ time_unit])
        #         t_col=ti/Tunits[time_unit]
        #     if ((vc - vm) * (vc0 - vm0) <= 0.0 ):
        #         #print("### RS:: time for vc=vm:", Rt[i-1], ti/Tunits[ time_unit])
        #         t_colmin=ti/Tunits[time_unit]
#-----------replace vm0 and vc0 for next loop
        vm0=vm
        vc0=vc
            

#        Rdec=(3.*Eiso/(Gm0**2* cgs.c**2 *4.*cgs.pi*n*cgs.mp))**(1/3.)



#-----------Pvm and Fvm
        Pvm=grb.Pvmax(Gm3,Bc)


#-----------Corrections on flux due to jet break effect
#            Fvm=grb.Fvmax(Pvm,n1,Ri,D28,zi)*fj*fDN /cgs.uJy

        for l,v in enumerate(vlist):

            Fvm=(1+zi)*Ne3*Pvm/(4.*cgs.pi*(D28*1e28)**2.) /cgs.uJy

#            Fvm_IC=1.
            fviewF=fview**3.

            if (t[i]<= txo):
                Fvti=fviewF *grb.Fv_sbpl(v/fview,va,vm,vc,Fvm,pp)  *Funit
            else:
                Fvti=fviewF *grb.Fv_RS(v/fview,va,vm,vc,Fvm,pp)  *Funit

#            if (v/fview>=vMax):
#                if (v/fview/vMax > 300.):
#                    Fvti=Fvti*np.exp(-300.+1.)
#                else:
#                    Fvti=Fvti*np.exp(-v/fview/vMax+1.)
#
            lgFv[l][i]=np.log10(Fvti)
#            lgFv[i]=np.log10(Fvti)
#            Fvt[i,0]=lgt[i]
#            Fvt[i,1]=lgFv[i]


#---------------Table for critical time
#         if (( t_dec <= to[i]/Tunits[ time_unit]) and ( t_dec > to[i-1]/Tunits[ time_unit])):
#             for j in range(0,inump):
#                  lgtdec[j]=np.log10( t_dec)
#                  lgFluxdec[j]=lgFv[i]-j/2.
# #                    lgFluxdec_n[j]=np.log10(fview_n**3. *Fvm_n)-j/2.
                
            
# #---------------Spectrum----------------------------
# #            mytime=Time[ Time_mode]
# #            if( (to[i]/cgs.day) >= mytime) and ( (to[i-1]/cgs.day) < mytime):
#         if( (t[i]/cgs.day) >= mytime) and ( (t[i-1]/cgs.day) < mytime):
#              nu_a=va/Nuunits[ nu_unit] *fview
#              nu_min=vm/Nuunits[ nu_unit] *fview
#              nu_cool=vc/Nuunits[ nu_unit] *fview

# #---------------Spectral index
#             dv = 1.05
#             Fvtip = fviewF *grb.Fv_sbpl(dv* v/fview,va,vm,vc,Fvm,pp)  *Funit
#             dF = Fvtip/Fvti
#             lgdF = math.log(dF,10)
#             lgdv = math.log(dv,10)
#             print("###[Spec Index]:",lgdF/lgdv)

# #---------------Table for critical syncrotron frequecies
#             Fvva=fviewF *grb.Fv_sbpl(va,va,vm,vc,Fvm,pp)
#             Fvvm=fviewF *grb.Fv_sbpl(vm,va,vm,vc,Fvm,pp)
#             Fvvc=fviewF *grb.Fv_sbpl(vc,va,vm,vc,Fvm,pp)
#             for j in range(0,inump):
#                  lgNua[j]=np.log10( nu_a)
#                  lgNum[j]=np.log10( nu_min)
#                  lgNuc[j]=np.log10( nu_cool)
#                  lgFluxa[j]=np.log10(Fvva)-j/2.
#                  lgFluxm[j]=np.log10(Fvvm)-j/2.
#                  lgFluxc[j]=np.log10(Fvvc)-j/2.

#             for j in range(1, inumv):
#                 vj=lgNu[j]
#                 Fvtj=fviewF *grb.Fv_sbpl(vj/fview,va,vm,vc,Fvm,pp)
#                 if (vj/fview>=vMax):
#                     if (vj/fview/vMax > 300.):
#                         Fvtj=Fvtj*np.exp(-300.+1.)
#                     else:
#                         Fvtj=Fvtj*np.exp(-vj/fview/vMax+1.)
#                 lgNu[j]=np.log10(vj/Nuunits[ nu_unit])
#                 lgFnu[j]=np.log10(Fvtj)





#---------------Parameters-----------
        para[i,0]=ti/Tunits[ time_unit]
        para[i,1]=Ri
        para[i,2]=Gm3
        para[i,3]=beta
        para[i,4]=n3
        para[i,5]=ei
        para[i,6]=Bc
        para[i,7]=Uph
        para[i,8]=UB
        para[i,9]=gm_a
        para[i,10]=gm_m
        para[i,11]=gm_c
        para[i,12]=gm_Max
        para[i,13]=va*fview
        para[i,14]=vm*fview
        para[i,15]=vc*fview
        para[i,16]=vMax*fview
        para[i,17]=Fvm
#            para[i,18]=Fvm_IC


#-----------small correction
    Gmt[0]=Gmt[1]
    lgt[0]=lgt[2]
    lgFv[:,0]=lgFv[:,2]
#    lgFv[0]=lgFv[2]
    # lgNu[0]=lgNu[2]
    # lgFnu[0]=lgFnu[2]
    lgt[1]=lgt[2]
    lgFv[:,1]=lgFv[:,2]
    # lgNu[1]=lgNu[2]
    # lgFnu[1]=lgFnu[2]
    for j in range(0,19):
        para[0,j]=para[1,j]

#-------transfer data to global
    paralist=para
    #  lgLC=np.column_stack((lgt,lgFv))
    #  lgSpec=np.column_stack((lgNu,lgFnu))
#         lgLC_SSC=np.column_stack((lgt,lgFv_SSC))
#         lgSpec_SSC=np.column_stack((lgNu,lgFnu_SSC))

#        np.savetxt( file_dir+'Output/vat_n.txt',vat) #t, va,vm,vc
#        np.savetxt( file_dir+'Output/Coolingt_n.txt',Coolingt)
#        np.savetxt( file_dir+'Output/Bt_n.txt',Bt_n)
#        np.savetxt( file_dir+'Output/para_n.txt',para) #t, Fvm, Gamma, R, Bc, gamma_m, gamma_c
#    print("### RS:: txo1,gm_cx1,gm_mx0,Rx,n1x0,n3x,e3x,Bcx,fn41x,Gm341x:", txo,gm_cx,gm_mx,Rx,n1x,n3x,e3x,Bcx,fn41x,Gm341x)
#        np.savetxt('/gtap/data/t.txt',10**(lgt))

#================end of Reverse shock==================
    lgFv = lgFv-3.0   # convert unit to "mJy"


#--------Estimate the flux of obs_time (add by fsy)
    estimate_model_log_flux = np.zeros(len(log_obs_time))
    for i,t in enumerate(log_obs_time):
        for iline in range(len(vlist)):
            if nu_obs[i] == vlist[iline]:
                vline = iline
                break
        model_log_time = lgt
        estimate_model_log_flux[i] = np.interp(t,model_log_time,lgFv[vline],left=-99.0,right=99.0)
    
    return 10**estimate_model_log_flux

