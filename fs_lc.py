# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 21:28:54 2017

@author: leiwh
"""
#-----------Scientific Library-------
import math
import numpy as np
import scipy.optimize

#-----------My Library---------
#Credit: Weihua Lei
import cgs
import astro
import grb



def FS_flux(z,E52,Gm0,n18,epsilon_B,epsilon_e,p, Band_mode, tobs):
#==============Forward shock=====================
#-------Input parameters: independent of model
        zi=z
        D28=astro.LUD(zi)/1.e28

        k=1.8

        n_dj=0.
        t_djs=0.
        t_dje=0.

        Lj50=0.
        q=0.1
        t_ej=1000.
        t_ejp=2000.
        t_ejf=5000.

        mytime=1.*cgs.day



        Band={'X-ray': 1e18,'Optical': 1e15,'Radio': 225e9, 'GeV': 2.4e23}
#        Time={'input':tday, 'adjust':10.0**lgtday}
        Dynmodel={'analytical': 1,'differential': 2}
        with_SSC={'Yes':1, 'No':0}
        with_jetbreak={'Yes':1, 'No':0}
        Tunits = {'Second': 1.0, 'Day': cgs.day}
        Nuunits= {'Hz': 1.0, 'keV': cgs.keV2Hz}

        dynamic_model='differential'
#        Band_mode='X-ray'
        v=Band[Band_mode]
        Funit=1.
        time_unit='Second'
        nu_unit='Hz'
        jet_break='Yes'
        inj_model='fallback'

        SSC = 'No'
        XIC = 'No'
        Smooth = 'No'


#-------refresh output parameters for jet
        t_dec=1.0
        t_Sed=1.0
        t_jet=1.0
        t_col=1.0
        t_min=1.0
        t_colmin=1.0
        nu_a=1.0
        nu_min=1.0
        nu_cool=1.0

#-------Input for jet

        n1=n18
        pp=p
        inum=500
        inumv=500
        inump=20

#-------Initialize the time, gamma, Fv arries in Log space
        t=np.logspace(-4., 3.0,inum)
        to=np.logspace(-4., 3.0,inum)  # t_observe

        Rt=np.logspace(14.,21.,inum)
        Gmt=np.logspace(2.,-5.,inum)


        para=np.zeros((inum,19))
        lgt=np.logspace(-4., 3.0,inum)
        lgFv=np.logspace(1.,6.,inum)
        lgFv_SSC=np.logspace(1.,6.,inum)

        lgNu=np.logspace(7.,30.,inumv)
        lgFnu=np.logspace(1.,6.,inumv)
        lgFnu_SSC=np.logspace(1.,6.,inumv)

#-------Initialize critical parameters
        lgNua=np.logspace(-4.,3.0,inump)
        lgNum=np.logspace(-4.,3.0,inump)
        lgNuc=np.logspace(-4.,3.0,inump)
        lgFluxa=np.logspace(-4.,3.0,inump)
        lgFluxm=np.logspace(-4.,3.0,inump)
        lgFluxc=np.logspace(-4.,3.0,inump)

        lgtdec=np.logspace(-4.,3.0,inump)
        lgFluxdec=np.logspace(-4.,3.0,inump)
        lgtSed=np.logspace(-4.,3.0,inump)
        lgFluxSed=np.logspace(-4.,3.0,inump)

        lgL=np.logspace(1.,5.,len(tobs))


#        lgLC=np.zeros((inum,2))

#-------start calculation
        Gmt[1]=Gm0

        R0=Rt[1]
        Eiso=1.e52 *E52
        t[0]=0.
        t[1]=R0/(2.*Gm0**2.0*cgs.c*cgs.day)
        to[0]=0.
        to[1]=t[1]

#        if (lgtday<math.log(t[1],10)):
#            lgtday=math.log(t[1],10)

#        Rdec=(3.*Eiso/(Gm0**2* cgs.c**2 *4.*cgs.pi*n*cgs.mp))**(1/3.)


#-------jet, LoS, view angle
        theta_obs_j=0.
        theta_j=4.445
        theta_obs=theta_obs_j
        thetaj=theta_j/cgs.deg
        thetaobs=theta_obs/cgs.deg
        fb=(1.-math.cos(thetaj))*2.*cgs.pi/(4.*cgs.pi)
        if (theta_obs <= theta_j):
            theta_view=0.
        else:
            theta_view=(thetaobs - thetaj)
        m0=1.e-5
        dm=0.0
        Ne0=0.0
        dNe=0.0
#-------true jet energy: E_j=E_iso*f_b
        Ej=Eiso*fb
        Mej0=Ej/((Gm0-1.)*cgs.c**2.0)
        epsilon_rad=0.0
        beta=1.0

        vm0=0.0
        vc0=0.0
        fview=1.
        fj=1.
        dE=0.0
#        m0=2.*cgs.pi*R0**3.*(1.-math.cos(theta_j))*n*cgs.mp/3.
        for i in range(1,inum):
            ni=n18*(Rt[i]/1.e18)**(-k)
            n1=ni
            if (i==1):
                Gmt[i]=Gm0
                beta=math.sqrt(1.-1/Gmt[i]**2.0)
                fview=(1.-beta)/(1.-beta*math.cos(theta_view) )
                t[i]=R0/(2.*Gm0**2.0*cgs.c) *(1.+zi)
                mt=m0
            else:
                dR=Rt[i]-Rt[i-1]
                dm=2.*cgs.pi*Rt[i]**2.*(1.-math.cos(thetaj))*n1*cgs.mp*dR
                mt=mt+dm
                if (Dynmodel[dynamic_model]==1):
#-------------------if (dynamic_model=="analytical")
                    Gmt[i] = grb.Gammat(Gm0,n1, Eiso, Rt[i],R0)
                else:
#-------------------if (dynamic_model=="differential")
                    tj=t[i-1]
                    if (inj_model == 'No'):
                        dE=0.
                    
                    elif ((tj > t_ej) and (tj< t_ejf)):
                        ts=t_ej
                        if (inj_model == 'fallback'):
                            tp=t_ejp
                            s=0.5
                            Lj=fb*Lj50*1.e50 *(0.5*((tj-ts)/(tp-ts))**(-0.5*s) + 0.5*((tj-ts)/(tp-ts))**(q*s) )**(-1./s)
                        else:
                            Lj=fb*Lj50*1.e50* (tj/ts)**(-q)

                        dE=(1.-beta)/beta/cgs.c/cgs.c/cgs.c * Lj*dR

                    dGm=-((Gmt[i-1]**2.-1.)*dm - dE)/(Mej0+mt*epsilon_rad+2.*(1.-epsilon_rad)*Gmt[i-1]*mt)
                    Gmt[i]=Gmt[i-1]+dGm
                beta=math.sqrt(1.-1/Gmt[i]**2.0)
                dtb=dR/(beta*cgs.c)
                dt=dtb/Gmt[i]/(Gmt[i]+math.sqrt(Gmt[i]**2.0-1.0)) *(1.+zi)
#---------------effect of on-axis and off-axis
                fview=(1.-beta)/(1.-beta*math.cos(theta_view) )
                dt=dt
                t[i]=t[i-1]+dt  #central engine time
            GM21=Gmt[i]
            ti=t[i]/fview  #off-viewer's time
            to[i]=ti
            Ri=Rt[i]
            lgt[i]=math.log(ti/Tunits[time_unit],10)

#-----------add density-jump
            if (n_dj>0.) and (ti>=t_djs) and (tj<=t_dje):
                n1=n_dj
            else:
                n1=ni


#-----------t_dec and t_Sed-------------------
            if ((mt-dm)*Gmt[i-1] <= Mej0) and (mt*Gmt[i] > Mej0):
#                print("### Deceleration:", Rt[i-1], grb.Rdec(Gm0,n1,Eiso), ti/Tunits[time_unit])
                t_dec=ti/Tunits[time_unit]

            if ((mt-dm) *cgs.c**2 <= Eiso*fb) and (mt*cgs.c**2 > Eiso*fb):
                Rsed0=grb.Rsed(n1,Eiso*fb)
#                if (Rt[i-1] <= Rsed0) and (Rt[i] > Rsed0):
#                print("### Sedov:", Rt[i-1], Rsed0, ti/Tunits[time_unit], Gmt[i])
                t_Sed=ti/Tunits[time_unit]

#-----------Internal energy---------
            ei=grb.e2(GM21,n1)
            n2=grb.n2_c(GM21,n1)

#-----------B field at comoving frame-------
            Bc=grb.Bco2(epsilon_B,ei)

#-----------Max Lorentz factor-------
            gm_Max = grb.gamma_Max(Bc)
#-----------Minimum Lorentz factor-------
            gm_m=grb.gamma_m2(epsilon_e,GM21,pp,gm_Max)
            gm_m0=gm_m
#-----------Cooling Lorentez factor
            gm_c=grb.gamma_c(Bc,GM21,ti*fview,zi)
            gm_c0=gm_c

#----------For deep Newtonian Phase----
#Introduce a factor fDN due to deep Newtonian
#See Sironi & Giannios, 2013, ApJ, 778, 107
            gm_syn = 1.
            fDN = 1.
            if (gm_m0 <= gm_c0):
                if (gm_m0 < gm_syn):
                    gm_m=gm_syn
                    fDN = gm_m0/gm_syn
                    if (pp > 3.):
                        fDN = (gm_m0/gm_syn)**((-1.+pp)/2.)
            else:
                if (gm_c0 < gm_syn):
                    gm_c0=gm_syn
                    fDN=gm_c0
#                    fDN_n = ((gm_syn_n -1.)/(gm_m_c -1.))**(1.-2.)


#-----------Y-parameter for SSC without K-N correction
            UB = Bc**2. /(8.*cgs.pi)
#            Coolingt[i,0]=ti/Tunits[time_unit]
#            Coolingt[i,2]=UB

            Y_SSC = 0.
            Uph=1.
            if SSC == 'Yes':
                if gm_c <= gm_m:
                    eta_e = 1.
                else:
                    eta_e = (gm_c/gm_m)**(2.-pp)
                ratio_epsilon = epsilon_e/epsilon_B
                Y_SSC = (-1. + math.sqrt(1.+4.* eta_e* ratio_epsilon))/2.
#                Coolingt[i,4]=UB* Y_SSC
                Uph=UB* Y_SSC
#                Coolingt[i,5]=Y_SSC

#-----------IC cooling of electrons by X-ray photons from internal shock, see Kumar et al. 2013
            if (XIC == 'Yes'):
#---------------X-ray flux from interal shock
                #Lx = 5.e46  # X-ray flux for SwJ1644+57
                #Lx = 10**48.2318 * t_n[i]**(-1.4)
                #Lx = 10**48.5 * (t_n[i]/cgs.day)**(-1.4)
                td_ce=t[i]/cgs.day #central engine time in unit of day
                if td_ce < 20.:
                    #if t_n[i]/cgs.day < 5.:
                        #Lx = 50.* 4. *10**46
                    #else:
                    Lx = 4. *10**46
                    #print t_n[i]/cgs.day, ', R= ',Ri_n, ', B=', Bc_n, 'n=', n1_n
                else:
                    Lx = 4. *10**46 * (td_ce /20.)**(-5/3.)
#---------------Comoving frame magnetic energy density and IC photons from interal shock
                Ux = Lx/(16.*cgs.pi* GM21**2. * Ri**2. *cgs.c)
#                UB_n = Bc_n**2. /(8.*cgs.pi)
                Y_XIC = Ux/UB
#                Coolingt_n[i,0]=t_n[i]/Tunits[time_unit]
#                Coolingt[i,1]=Ux
#                Coolingt_n[i,2]=UB_n
#                Coolingt[i,3]=Y_XIC
                #gm_c_n=grb.gamma_c(Bc_n,GM21_n,ti_n,zi) /(1. + Y_n)

#---------------Include the coorections due to the K-N cross-section
                betax = 0.7
                nux = 10.* 1.e-3 *(1.+zi)/GM21
                rex = (cgs.mecc2MeV/nux)**(1.-betax)
                Cx = Y_XIC* rex

                def  fkn(ge): return Cx* ge**betax + ge* (1.+ Y_SSC) -gm_c

                if Cx >= (gm_c -1.):
                    gm_c_kn =1.
                else:
                #gm_c_n_kn = scipy.optimize.newton(fkn,1.)
                    gm_c_kn = scipy.optimize.brenth(fkn,1., gm_c)
                #gm_c_n_kn = gm_c_n/(1.+ Y_n* rex)
                #print gm_c_n_kn
                gm_c = gm_c_kn
            else:
                gm_c = gm_c/(1.+ Y_SSC)

#-----------v_Max, v_m, v_c, v_a----------
            vMax = grb.nu_gme(GM21,gm_Max, Bc,zi)
            vm=grb.nu_gme(GM21,gm_m, Bc,zi)
            vc=grb.nu_gme(GM21,gm_c, Bc,zi)
            va=grb.nu_aB(Bc,pp,n1,vm,vc,GM21,Ri,zi)
            gm_a=(va*(2.*cgs.pi *cgs.me *cgs.c)*(1.+zi)/(GM21 *cgs.qe *Bc))**0.5

#            print 't=',t_n[i]/cgs.day, ', vMax=',vMax_n,', vm=', vm_n, ', vc=',vc_n

#            vat[i,0]=np.log10(ti/Tunits[time_unit])
#            vat[i,1]=np.log10(va*fview)
#            vat[i,2]=np.log10(vm*fview)
#            vat[i,3]=np.log10(vc*fview)
#            vat[i,4]=np.log10(vMax*fview)

#-----------t_min, t_col and t_colmin, t_jet
            if (i>1):
                if ((v - vm) * (v - vm0) <= 0.0 ):
#                    print("### time for vm cross:", Rt[i-1], ti/Tunits[time_unit])
                    t_min=ti/Tunits[time_unit]
                if ((v - vc) * (v - vc0) <= 0.0 ):
#                    print("### time for vc cross:", Rt[i-1], ti/Tunits[time_unit])
                    t_col=ti/Tunits[time_unit]
                if ((vc - vm) * (vc0 - vm0) <= 0.0 ):
#                    print("### time for vc=vm:", Rt[i-1], ti/Tunits[time_unit])
                    t_colmin=ti/Tunits[time_unit]
#-----------replace vm0 and vc0 for next loop
            vm0=vm
            vc0=vc
#-----------Jet break time: t_jet
            if (with_jetbreak[jet_break]==1):
                if (thetaj < (1./GM21) ):
                    if (fj == 1.):
#                        print("### jet break time:", ti/Tunits[time_unit])
                        t_jet=ti/Tunits[time_unit]
                    fj=thetaj**2. /(1./GM21**2.)

#-----------Corrections on off-axis jet for flux, see Salafia et al. 2016, arXiv:1601.03735
            if (thetaj<1./GM21):
                thetaj = 1./GM21

            thetajs=thetaj-1/GM21
            if (thetaobs <= thetajs):
                fviewF=1.
            elif ((thetaobs > thetajs) and (thetaobs <= thetaj)):
                fviewF = 1. - GM21*(thetaobs-thetajs)/2.
            else:
                fviewF=0.5* fview**(3.-np.sqrt(2.)*thetaj**(1/3.))


#-----------Pvm and Fvm
            Pvm=grb.Pvmax(GM21,Bc)

#-----------Corrections on flux due to jet break effect
            Ne=mt/(cgs.mp*fb) 
            Fvm=(1+zi)*Ne*Pvm/(4.*cgs.pi*(D28*1e28)**2.)*fj*fDN /cgs.uJy

#            Fvm=grb.Fvmax(Pvm,n1,Ri,D28,zi)*fj*fDN /cgs.uJy

            Fvm_IC=1.
            fviewF=fview**3.

#            Fvti=fviewF *grb.Fv(v/fview,va,vm,vc,Fvm,pp)  *Funit
            if (Smooth=='Yes'):
                Fvti=fviewF *grb.Fv(v/fview,va,vm,vc,Fvm,pp)  *Funit
            else:
                Fvti=fviewF *grb.Fv_sbpl(v/fview,va,vm,vc,Fvm,pp)  *Funit

            if (v/fview>=vMax):
                if (v/fview/vMax > 300.):
                    Fvti=Fvti*np.exp(-300.+1.)
                else:
                    Fvti=Fvti*np.exp(-v/fview/vMax+1.)

            lgFv[i]=np.log10(Fvti)


            if (with_SSC[SSC]==1):
                Fvm_IC = grb.Fvmax_IC(Fvm, n1, Ri)
#                Fvt_IC_n = fview_n**3. *grb.FvIC_simp(v/fview_n,gm_m_n, gm_c_n, vm_n, vc_n, Fvm_IC_n, pp_n)  *Funit
#                Fvt_IC = fview**3. *grb.FvIC(v/fview,gm_a,gm_m, gm_c, va,vm, vc, vMax, Fvm_IC, GM21,pp) *Funit
                Fvt_IC = fviewF *grb.FvIC(v/fview,gm_a,gm_m, gm_c, va,vm, vc, vMax, Fvm_IC, GM21,pp) *Funit
                lgFv_SSC[i] = math.log(Fvt_IC,10)
#                LC[i,2]=lgFv_SSC[i]


#---------------Table for critical time
            if ((t_dec <= to[i]/Tunits[time_unit]) and (t_dec > to[i-1]/Tunits[time_unit])):
                for j in range(0,inump):
                    lgtdec[j]=np.log10(t_dec)
                    lgFluxdec[j]=lgFv[i]-j/2.
#                    lgFluxdec_n[j]=np.log10(fview_n**3. *Fvm_n)-j/2.
            if ((t_Sed <= to[i]/Tunits[time_unit]) and (t_Sed > to[i-1]/Tunits[time_unit])):
                for j in range(0,inump):
                    lgtSed[j]=np.log10(t_Sed)
                    lgFluxSed[j]=lgFv[i]-j/2.

#---------------Spectrum----------------------------
#            mytime=Time[Time_mode]
            if( (to[i]/cgs.day) >= mytime) and ( (to[i-1]/cgs.day) < mytime):
                nu_a=va/Nuunits[nu_unit] *fview
                nu_min=vm/Nuunits[nu_unit] *fview
                nu_cool=vc/Nuunits[nu_unit] *fview
#---------------Spectral index
#                dv = 1.05
#                if (self.Smooth=='Yes'):
#                    Fvtip=fviewF *grb.Fv(dv* v/fview,va,vm,vc,Fvm,pp) *Funit
#                else:
#                    Fvtip=fviewF *grb.Fv_sbpl(dv* v/fview,va,vm,vc,Fvm,pp) *Funit
#                dF = Fvtip/Fvti
#                lgdF = math.log(dF,10)
#                lgdv = math.log(dv,10)
#                print "###[Spec Index]:",lgdF/lgdv

#---------------Table for critical syncrotron frequecies
                if (Smooth=='Yes'):
                    Fvva=fviewF *grb.Fv(va,va,vm,vc,Fvm,pp)
                    Fvvm=fviewF *grb.Fv(vm,va,vm,vc,Fvm,pp)
                    Fvvc=fviewF *grb.Fv(vc,va,vm,vc,Fvm,pp)
                else:
                    Fvva=fviewF *grb.Fv_sbpl(va,va,vm,vc,Fvm,pp)
                    Fvvm=fviewF *grb.Fv_sbpl(vm,va,vm,vc,Fvm,pp)
                    Fvvc=fviewF *grb.Fv_sbpl(vc,va,vm,vc,Fvm,pp)
                
                for j in range(0,inump):
                    lgNua[j]=np.log10(nu_a)
                    lgNum[j]=np.log10(nu_min)
                    lgNuc[j]=np.log10(nu_cool)
                    lgFluxa[j]=np.log10(Fvva)-j/2.
                    lgFluxm[j]=np.log10(Fvvm)-j/2.
                    lgFluxc[j]=np.log10(Fvvc)-j/2.

                for j in range(1, inumv):
                    vj=lgNu[j]
                    if (Smooth=='Yes'):
                        Fvtj=fviewF *grb.Fv(vj/fview,va,vm,vc,Fvm,pp)
                    else:
                        Fvtj=fviewF *grb.Fv_sbpl(vj/fview,va,vm,vc,Fvm,pp)

                    if (vj/fview>=vMax):
                        if (vj/fview/vMax > 300.):
                            Fvtj=Fvtj*np.exp(-300.+1.)
                        else:
                            Fvtj=Fvtj*np.exp(-vj/fview/vMax+1.)
#                    vmax=grb.nu_Max(GM21,zi)
#                    Fvtj=grb.Fv_simp(vj,vm,vc,vmax, Fvm,pp)
                    lgNu[j]=np.log10(vj/Nuunits[nu_unit])
                    lgFnu[j]=np.log10(Fvtj)
#                    lgvFluxv[j]=np.log10(vj*Fvtj)
#                    Fvv[j,0]=lgNuv[j]
#                    Fvv[j,1]=lgFluxv[j]
#                    Spec[j,0]=lgNuv[j]
#                    Spec[j,1]=lgFluxv[j]
#                    Spec[j,2]=lgvFluxv[j]

                    if (with_SSC[SSC]==1):
#                        Fvtj_IC= fview_n**3. *grb.FvIC_simp(vj_n/fview_n,gm_m_n, gm_c_n, vm_n, vc_n, Fvm_IC_n, pp_n)
#                        Fvtj_IC =fview**3. *grb.FvIC(vj/fview,gm_a, gm_m, gm_c,va,vm, vc, vMax, Fvm_IC, GM21, pp)
                        Fvtj_IC =fviewF *grb.FvIC(vj/fview,gm_a, gm_m, gm_c,va,vm, vc, vMax, Fvm_IC, GM21, pp)
                        lgFnu_SSC[j]=np.log10(Fvtj_IC)
#                        lgvFluxv_SSC[j]=np.log10(vj*Fvtj_IC)
#                        Spec[j,3]=lgFluxv_SSC[j]
#                        Spec[j,4]=lgvFluxv_SSC[j]


#---------------Parameters-----------
            para[i,0]=ti/Tunits[time_unit]
            para[i,1]=Ri
            para[i,2]=GM21
            para[i,3]=beta
            para[i,4]=n2
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
            para[i,18]=Fvm_IC


#-----------small correction
        Gmt[0]=Gmt[1]
        lgt[0]=lgt[2]
        lgFv[0]=lgFv[2]
        lgNu[0]=lgNu[2]
        lgFnu[0]=lgFnu[2]
        lgt[1]=lgt[2]
        lgFv[1]=lgFv[2]
        lgNu[1]=lgNu[2]
        lgFnu[1]=lgFnu[2]
        lgFv_SSC[0]=lgFv_SSC[2]
        lgFnu_SSC[0]=lgFnu_SSC[2]
        lgFv_SSC[1]=lgFv_SSC[2]
        lgFnu_SSC[1]=lgFnu_SSC[2]
        for j in range(0,19):
            para[0,j]=para[1,j]

#-------transfer data to global
        paralist=para
        lgLC=np.column_stack((lgt,lgFv))
        lgSpec=np.column_stack((lgNu,lgFnu))
        lgLC_SSC=np.column_stack((lgt,lgFv_SSC))
        lgSpec_SSC=np.column_stack((lgNu,lgFnu_SSC))


        for i in range(0,len(tobs)):
             listt=lgt
             listt=listt.tolist()

             to=tobs[i]
#             fn=lgFv_n[i]
             listt.append(to)
             listt.sort()
             indexo=listt.index(to)
             if (indexo == 0):
                 j=indexo+1
             else:
                 j=indexo
                 if (indexo==len(lgt)):
                        j=indexo-1
#-------------------Estimate Fv_w at time t_n
             lgL[i]=(to-lgt[j-1])*(lgFv[j]-lgFv[j-1])/(lgt[j]-lgt[j-1]) + lgFv[j-1]

        lgL=lgL+np.log10(v*cgs.uJy)

        return lgL


#        return lgLC

#        np.savetxt(file_dir+'Output/vat_n.txt',vat) #t, va,vm,vc
#        np.savetxt(file_dir+'Output/Coolingt_n.txt',Coolingt)
#        np.savetxt(file_dir+'Output/Bt_n.txt',Bt_n)
#        np.savetxt(file_dir+'Output/para_n.txt',para) #t, Fvm, Gamma, R, Bc, gamma_m, gamma_c

#================end of Forward shock==================
