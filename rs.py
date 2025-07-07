
#================================================
#==============Reverse shock=====================

#-----------Scientific Library-------
import math
import numpy as np
import scipy.optimize

#-----------My Library---------
#Credit: Weihua Lei
try:
    from . import cgs,astro,grb
except ImportError:
    import cgs,astro,grb


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
    RS_Correction='No'  # Correct for the RS
    xi_s=1/3.**0.5  #crrection for radial spreading of RS
    xi_gm34=7/4.0 # correction for gamma34 in RS

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

#-------jet, LoS, view angle
    #theta_obs=theta_obs
    thetaj=theta_j/cgs.deg
    thetaobs=theta_obs/cgs.deg
    fb=(1.-math.cos(thetaj))*2.*cgs.pi/(4.*cgs.pi)

    Gm4=Gm0
    Gm3x=Gm0
    Gm3=Gm0-1.e-4
    #Delta0=t90*cgs.c/(1.+zi)

    An0=n18 *1.e18**k
    Ne0=Eiso/(Gm0*cgs.mp*cgs.c**2.)
    xi=0.

    Delta0=R0/(3**0.5 * Gm0**2.)
    Delta=Delta0

    #spreading radius
    Rs=3**0.5 * Gm0**2. *Delta0
    beta4=math.sqrt(1.-1/Gm0**2.0)


    fview=1.
    Deltax=Rt[1]/(3**0.5 * Gm0**2.)  

    for i in range(1,inum):
        ni= n18*(Rt[i]/1.e18)**(- k)
        n1=ni

        Ri=Rt[i]
        dR=Ri-Rt[i-1]
        t[i]=Ri/(2.*Gm3**2.0 *cgs.c) *(1.+zi)
        ti=t[i]

        # check shell spreading
#        if (Ri<Rs):
#            Delta=Delta0
#        else:
#            Delta=Ri/(3**0.5 * Gm0**2.)

        Delta=Ri/(3**0.5 * Gm0**2.)   

        # check RS crossing
        if (xi< Deltax): #before crossing
            n4=Eiso/(Gm0*cgs.mp*cgs.c**2. *4.*cgs.pi*Ri**2. *Gm0*Delta)
            fn41=n4/n1
            beta3=math.sqrt(1.-1/Gm3**2.0)
            #Gm34=Gm0*Gm3*(1.-beta3*beta4)
            #Gm34=(Gm0/Gm3+Gm3/Gm0)/2.
            #Gm34=0.5*Gm3**2./fn41+ 1.
            Gm34=((Gm3**2.-1.)/fn41+ 1.)**0.5
            #print("Gm3=,Gm34=",Gm3,Gm34)
            gmcat3 = (4.*Gm34+1.)/(3.*Gm34)
            n43=(gmcat3-1.)/(gmcat3*Gm34+1.)
            n3=n4/n43
#            dx=dR/(Gm0*fn41**0.5 * (1.-1/4.))
            dx=dR/(Gm0*fn41**0.5 * (1.-Gm0/Gm3 * n43))

            xi=xi+dx
            xp=Gm0/Gm3 * n43* xi
            Ne3=n3*4.*cgs.pi*Ri**2.* Gm3*xp
            e3=(Gm34-1.)*n3*cgs.mp*cgs.c**2.
            Bc=(8.*cgs.pi *e3 * epsilon_B)**0.5

#-----------Max Lorentz factor-------
            gm_Max = grb.gamma_Max(Bc)
#-----------Minimum Lorentz factor-------
            gm_m=grb.gamma_m2(epsilon_e,Gm34,pp,gm_Max)
#-----------Cooling Lorentez factor
            gm_c=grb.gamma_c(Bc,Gm3,ti*fview,zi)

            #parameters at crossing
            e3x=e3
            n3x=n3
            n1x=n1
            fn41x=fn41
            Bcx=Bc
            txo=ti
            Rx=Ri
            Gm3x=Gm3
            Gm34x=Gm34
            gm_mx=gm_m
            gm_cx=gm_c
            Deltax=Delta

            #Gm3=(fn41*(Gm34**2.-1.) +1.)**0.5
            #Gm3=Gm0-1.e-4
            Gm3=grb.Gammat(Gm4,ni, Eiso, Ri,R0)

            #print("Gm34=, Gm3=, Delta=, Delta0=, xi=",Gm34,Gm3,Delta,Delta0,fn41**0.5/Gm0)

        else: #after crossing
            Ne3=Ne0

            gthick=(7-2.*k)/2.
            gthin=2-k/2.
            g= ((Gm34-1.)*gthick+gthin)/Gm34
            Gm3=Gm3x*(Ri/Rx)**(-g)
            n3=n3x*(Ri/Rx)**(-2.*g)

            if (Gm3 <=1.):
                Gm3=1.+1.e-6

            gthick_e=(26-4.*k)/3.
            gthin_e=4.-k
            g_e= ((Gm34-1.)*gthick_e+gthin_e)/Gm34
            e3=e3x*(Ri/Rx)**(-g_e)

            Bc=(8.*cgs.pi *e3 * epsilon_B)**0.5

            gm_Max = grb.gamma_Max(Bc)
            gm_m=gm_mx*(Ri/Rx)**(-g_e+2*g)
            gm_c=gm_cx*(Ri/Rx)**(-g_e+2*g)

            Pvmx=grb.Pvmax(Gm3x,Bcx)

            #print("### RS:: tx,gm_cx0,gm_mx0,Rx,n1x0,n3x,e3x,Bcx,fn41x,Gm341x:", txo,gm_cx,gm_mx,Rx,n1x,n3x,e3x,Bcx,fn41x,Gm34x)



        #print("Gmma3=",Gm3x)
        theta_view=grb.theta_off(Gm3,thetaobs,thetaj)
        fview=grb.aoff(Gm3,theta_view)


        to[i]=ti/fview  #off-viewer's time
        lgt[i]=math.log(ti/Tunits[ time_unit],10)
        ei=e3

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
#-----------Corrections on off-axis jet for flux, see Beniamini et al. 2023
            fviewF=grb.FvOff(Gm3,thetaobs,thetaj)
            #fviewF=fview**3.

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
        para[i,3]=beta3
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
        para[i,18]=Gm34


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

