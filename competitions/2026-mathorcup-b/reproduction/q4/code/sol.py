from pathlib import Path
import json, math, time
import numpy as np
from scipy.integrate import solve_ivp

G = 9.81
DRY = 1e-4
ORDER = 1
FAST = None


import ctypes,hashlib,subprocess,sys,json,re
from pathlib import Path
import numpy as np

SOURCE=r'''
#include <math.h>
#include <stdlib.h>
#include <string.h>
static const double G=9.81;
static double wallp(double h,double w){
  if(h<=0)return 0.;double c=sqrt(G*h);
  if(w<=0){double cs=fmax(0.,c+.5*w),hs=cs*cs/G;return .5*G*hs*hs;}
  double f=w/c,lo=1.,hi=1.+2.*f+f*f;
  for(int j=0;j<60;j++){double r=.5*(lo+hi),value=(r-1.)*sqrt((r+1.)/(2.*r));if(value<f)lo=r;else hi=r;}
  double hs=h*.5*(lo+hi);return .5*G*hs*hs;
}
static double mm(double a,double b,double c){return a*b>0&&a*c>0?copysign(fmin(fmin(fabs(a),fabs(b)),fabs(c)),a):0.;}
static void state(const double *w,const double *s,double sign,double *u,double *z){
  double h=fmax(0.,w[0]+sign*s[0]-(w[3]+sign*s[3]));
  *z=w[3]+sign*s[3];u[0]=h;u[1]=h*(w[1]+sign*s[1]);u[2]=h*(w[2]+sign*s[2]);
}
static void hll(const double *L,const double *R,double zl,double zr,int axis,double dry,double *F,double *CL,double *CR){
  double vl[2]={0,0},vr[2]={0,0},UL[3],UR[3],FL[3],FR[3];
  for(int k=0;k<2;k++){if(L[0]>dry)vl[k]=L[k+1]/L[0];if(R[0]>dry)vr[k]=R[k+1]/R[0];}
  double zs=fmax(zl,zr),hl=fmax(0.,L[0]+zl-zs),hr=fmax(0.,R[0]+zr-zs);
  double nl=vl[axis],nr=vr[axis],cl=sqrt(G*hl),cr=sqrt(G*hr);
  double sl=fmin(nl-cl,nr-cr),sr=fmax(nl+cl,nr+cr);
  if(hr<=0.){sl=nl-cl;sr=nl+2*cl;}if(hl<=0.){sl=nr-2*cr;sr=nr+cr;}
  UL[0]=hl;UR[0]=hr;for(int k=0;k<2;k++){UL[k+1]=hl*vl[k];UR[k+1]=hr*vr[k];}
  for(int k=0;k<3;k++){FL[k]=UL[k]*nl;FR[k]=UR[k]*nr;}
  FL[axis+1]+=.5*G*hl*hl;FR[axis+1]+=.5*G*hr*hr;
  for(int k=0;k<3;k++){double mid=sr-sl>dry?(sr*FL[k]-sl*FR[k]+sl*sr*(UR[k]-UL[k]))/(sr-sl):0.;F[k]=sl>=0?FL[k]:(sr<=0?FR[k]:mid);}
  *CL=.5*G*(L[0]*L[0]-hl*hl);*CR=.5*G*(R[0]*R[0]-hr*hr);
}
int fluxes(int ny,int nx,int order,int bc,double dry,const double *V,const double *Z,
 double *FX,double *FY,double *XL,double *XR,double *YL,double *YR,double *SX,double *SY,double *entryH,double *entryZ){
  int nyp=ny+2,nxp=nx+2,n=nyp*nxp;double *W=calloc(n*4,sizeof(double)),*sx=calloc(n*4,sizeof(double)),*sy=calloc(n*4,sizeof(double));
  if(!W||!sx||!sy){free(W);free(sx);free(sy);return 1;}
  for(int p=0;p<n;p++){double h=V[3*p];W[4*p]=h+Z[p];W[4*p+3]=Z[p];if(h>dry){W[4*p+1]=V[3*p+1]/h;W[4*p+2]=V[3*p+2]/h;}}
  if(order==2){
    for(int i=0;i<nyp;i++)for(int j=1;j<nxp-1;j++){
      int p=i*nxp+j;for(int k=0;k<4;k++)sx[4*p+k]=mm(2*(W[4*p+k]-W[4*(p-1)+k]),.5*(W[4*(p+1)+k]-W[4*(p-1)+k]),2*(W[4*(p+1)+k]-W[4*p+k]));
      double h=V[3*p];if(h<=dry||V[3*(p-1)]<=dry||V[3*(p+1)]<=dry||h<.5*fabs(sx[4*p]-sx[4*p+3]))for(int k=0;k<4;k++)sx[4*p+k]=0.;
    }
    for(int i=1;i<nyp-1;i++)for(int j=0;j<nxp;j++){
      int p=i*nxp+j;for(int k=0;k<4;k++)sy[4*p+k]=mm(2*(W[4*p+k]-W[4*(p-nxp)+k]),.5*(W[4*(p+nxp)+k]-W[4*(p-nxp)+k]),2*(W[4*(p+nxp)+k]-W[4*p+k]));
      double h=V[3*p];if(h<=dry||V[3*(p-nxp)]<=dry||V[3*(p+nxp)]<=dry||h<.5*fabs(sy[4*p]-sy[4*p+3]))for(int k=0;k<4;k++)sy[4*p+k]=0.;
    }
    if(bc==1||bc==2){
      for(int i=0;i<nyp;i++)for(int k=0;k<4;k++){
        sx[4*(i*nxp)+k]=sx[4*(i*nxp+(bc==1?1:nx))+k];
        sx[4*(i*nxp+nx+1)+k]=sx[4*(i*nxp+(bc==1?nx:1))+k];
      }
      for(int j=0;j<nxp;j++)for(int k=0;k<4;k++){
        sy[4*j+k]=sy[4*((bc==1?1:ny)*nxp+j)+k];
        sy[4*((ny+1)*nxp+j)+k]=sy[4*((bc==1?ny:1)*nxp+j)+k];
      }
    }
  }
  for(int i=0;i<ny;i++)for(int j=0;j<=nx;j++){
    int p=(i+1)*nxp+j,q=p+1,f=i*(nx+1)+j;double L[3],R[3],zl,zr;
    if(order==2){state(W+4*p,sx+4*p,.5,L,&zl);state(W+4*q,sx+4*q,-.5,R,&zr);}
    else{memcpy(L,V+3*p,3*sizeof(double));memcpy(R,V+3*q,3*sizeof(double));zl=Z[p];zr=Z[q];}
    if(j==0&&bc==3){memcpy(L,R,3*sizeof(double));L[1]*=-1;zl=zr;}
    if(j==nx&&(bc==3||bc==4)){memcpy(R,L,3*sizeof(double));R[1]*=-1;zr=zl;}
    hll(L,R,zl,zr,0,dry,FX+3*f,XL+f,XR+f);
    if(j==0&&bc==3){FX[3*f]=0;FX[3*f+2]=0;FX[3*f+1]=wallp(R[0],R[0]>dry?-R[1]/R[0]:0.);}
    if(j==nx&&(bc==3||bc==4)){FX[3*f]=0;FX[3*f+2]=0;FX[3*f+1]=wallp(L[0],L[0]>dry?L[1]/L[0]:0.);}
    if(j==nx){for(int k=0;k<3;k++)entryH[3*i+k]=L[k];entryZ[i]=zl;}
  }
  for(int i=0;i<=ny;i++)for(int j=0;j<nx;j++){
    int p=i*nxp+j+1,q=p+nxp,f=i*nx+j;double L[3],R[3],zl,zr;
    if(order==2){state(W+4*p,sy+4*p,.5,L,&zl);state(W+4*q,sy+4*q,-.5,R,&zr);}
    else{memcpy(L,V+3*p,3*sizeof(double));memcpy(R,V+3*q,3*sizeof(double));zl=Z[p];zr=Z[q];}
    if(i==0&&bc==3){memcpy(L,R,3*sizeof(double));L[2]*=-1;zl=zr;}
    if(i==ny&&bc==3){memcpy(R,L,3*sizeof(double));R[2]*=-1;zr=zl;}
    hll(L,R,zl,zr,1,dry,FY+3*f,YL+f,YR+f);
    if(i==0&&bc==3){FY[3*f]=0;FY[3*f+1]=0;FY[3*f+2]=wallp(R[0],R[0]>dry?-R[2]/R[0]:0.);}
    if(i==ny&&bc==3){FY[3*f]=0;FY[3*f+1]=0;FY[3*f+2]=wallp(L[0],L[0]>dry?L[2]/L[0]:0.);}
  }
  if(order==2)for(int i=0;i<ny;i++)for(int j=0;j<nx;j++){
    int p=(i+1)*nxp+j+1,f=i*nx+j;double L[3],R[3],zl,zr;
    state(W+4*p,sx+4*p,-.5,L,&zl);state(W+4*p,sx+4*p,.5,R,&zr);SX[f]=-.5*G*(L[0]+R[0])*(zr-zl);
    state(W+4*p,sy+4*p,-.5,L,&zl);state(W+4*p,sy+4*p,.5,R,&zr);SY[f]=-.5*G*(L[0]+R[0])*(zr-zl);
  }
  free(W);free(sx);free(sy);return 0;
}
'''


class FastFlux:
    def __init__(self,build_dir):
        self.sha=hashlib.sha256(SOURCE.encode()).hexdigest();out=Path(build_dir)/self.sha[:16];out.mkdir(parents=True,exist_ok=True)
        src=out/'flux.c';lib=out/'flux.so'

        if True:
            src.write_text(SOURCE)
            command=['cc','-O3','-std=c99','-shared','-fPIC',str(src),'-o',str(lib),'-lm']
            if sys.platform=='darwin':
                probe=subprocess.run(['xcrun','--show-sdk-path'],capture_output=True,text=True,check=True)
                default=Path(probe.stdout.strip());candidates=[]
                for candidate in default.parent.glob('MacOSX*.sdk'):
                    match=re.fullmatch(r'MacOSX(\d+)\.(\d+)\.sdk',candidate.name)
                    if match:candidates.append(((int(match[1]),int(match[2])),candidate))


                sdk=min(candidates,key=lambda x:x[0])[1] if candidates else default
                command+=['-isysroot',str(sdk)]
            run=subprocess.run(command,capture_output=True,text=True)
            (out/'build_command.json').write_text(json.dumps(command))
            (out/'build.stdout').write_text(run.stdout);(out/'build.stderr').write_text(run.stderr)
            if run.returncode:raise RuntimeError('C99 flux compilation failed: '+run.stderr)
        self.command=json.loads((out/'build_command.json').read_text())
        self.library_sha256=hashlib.sha256(lib.read_bytes()).hexdigest()
        self.library=ctypes.CDLL(str(lib.resolve()));self.fn=self.library.fluxes
        ptr=ctypes.POINTER(ctypes.c_double)
        self.fn.argtypes=[ctypes.c_int]*4+[ctypes.c_double]+[ptr]*12
        self.fn.restype=ctypes.c_int

    def fluxes(self,V,Z,order,boundary,east_wall=False):
        if np.iscomplexobj(V) or np.iscomplexobj(Z):raise TypeError('real-valued states required')
        if np.ma.isMaskedArray(V) or np.ma.isMaskedArray(Z):raise TypeError('masked arrays must be resolved before the native boundary')
        V=np.asarray(V,dtype=np.float64,order='C');Z=np.asarray(Z,dtype=np.float64,order='C')
        if V.ndim!=3 or V.shape[-1]!=3 or Z.ndim!=2 or V.shape[:2]!=Z.shape or min(Z.shape)<3:raise ValueError('padded V shape (ny+2,nx+2,3) and matching bed required')
        if not np.isfinite(V).all() or not np.isfinite(Z).all() or np.any(V[...,0]<0):raise ValueError('finite nonnegative-thickness states and finite bed required')
        if order not in (1,2) or boundary not in ('wall','open','outflow','slope','periodic'):raise ValueError('unsupported flux scheme/boundary')
        if not np.isfinite(self.dry) or self.dry<=0:raise ValueError('positive finite dry threshold required')
        if east_wall and boundary in ('slope','periodic'):raise ValueError('inlet combined with slope/periodic is not supported by the native kernel')
        ny,nx=Z.shape[0]-2,Z.shape[1]-2
        fx=np.zeros((ny,nx+1,3));fy=np.zeros((ny+1,nx,3));xl=np.zeros((ny,nx+1));xr=xl.copy();yl=np.zeros((ny+1,nx));yr=yl.copy()
        sx=np.zeros((ny,nx));sy=sx.copy();eh=np.zeros((ny,3));ez=np.zeros(ny)
        arrays=[np.ascontiguousarray(V),np.ascontiguousarray(Z),fx,fy,xl,xr,yl,yr,sx,sy,eh,ez]
        bc=3 if boundary=='wall' else (4 if east_wall else (1 if boundary=='slope' else (2 if boundary=='periodic' else 0)))
        status=self.fn(ny,nx,order,bc,self.dry,*[a.ctypes.data_as(ctypes.POINTER(ctypes.c_double)) for a in arrays])
        if status:raise MemoryError('C flux allocation failed')
        left=eh[:,None,:]
        return fx,fy,xl,xr,yl,yr,sx,sy,left,ez[:,None]

import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.optimize import minimize


def matern52(x,y,amplitude,length):
    r=np.sqrt(np.sum(((np.asarray(x)[:,None,:]-np.asarray(y)[None,:,:])/length)**2,axis=-1))
    a=np.sqrt(5)*r
    return amplitude**2*(1+a+a*a/3)*np.exp(-a)


def morris_design(dimension=4,trajectories=8,levels=6,seed=83031):
    rng=np.random.default_rng(seed);delta=levels/(2*(levels-1));points=[];orders=[]
    for _ in range(trajectories):
        direction=rng.choice([-1,1],dimension)
        base=rng.integers(0,levels//2,size=dimension)/(levels-1)
        x=np.where(direction>0,base,1-base)
        order=rng.permutation(dimension);track=[x.copy()]
        for j in order:
            x=x.copy();x[j]+=delta*direction[j];track.append(x.copy())
        points.append(track);orders.append(order)
    return np.array(points),np.array(orders)


def morris_effects(points,responses):
    points=np.asarray(points);responses=np.asarray(responses);ee=np.full((len(points),points.shape[2]),np.nan)
    for r,trajectory in enumerate(points):
        for j in range(len(trajectory)-1):
            difference=trajectory[j+1]-trajectory[j];changed=np.flatnonzero(abs(difference)>1e-10)
            if len(changed)!=1:raise ValueError('Morris step must change exactly one coordinate')
            d=changed[0];ee[r,d]=(responses[r,j+1]-responses[r,j])/difference[d]
    return ee,{'mu':np.mean(ee,axis=0).tolist(),'mu_star':np.mean(abs(ee),axis=0).tolist(),'sigma':np.std(ee,axis=0,ddof=1).tolist()}


class ARGP:
    def fit(self,xlow,ylow,xhigh,yhigh):
        ylow=np.asarray(ylow);yhigh=np.asarray(yhigh)
        self.xl=np.asarray(xlow);self.xh=np.asarray(xhigh)
        if self.xl.ndim!=2 or self.xh.ndim!=2 or self.xl.shape[1]!=self.xh.shape[1] or len(ylow)!=len(self.xl) or len(yhigh)!=len(self.xh) or not all(np.isfinite(v).all() for v in [self.xl,self.xh,ylow,yhigh]):raise ValueError('finite matching low/high GP data required')
        self.ml=float(np.mean(ylow));self.mh=float(np.mean(yhigh))
        self.scale=max(float(np.std(np.r_[ylow,yhigh])),1e-3)
        self.y=np.r_[(ylow-self.ml)/self.scale,(yhigh-self.mh)/self.scale]
        self.nL=len(xlow);self.dim=self.xl.shape[1];self.noise=1e-8
        def unpack(theta):
            k=self.dim
            return np.exp(theta[0]),np.exp(theta[1:k+1]),np.exp(theta[k+1]),np.exp(theta[k+2:2*k+2]),theta[-1]
        self.unpack=unpack
        def matrix(theta):
            al,ll,ad,ld,rho=unpack(theta)
            Lll=matern52(self.xl,self.xl,al,ll)
            Llh=rho*matern52(self.xl,self.xh,al,ll)
            Hhh=rho*rho*matern52(self.xh,self.xh,al,ll)+matern52(self.xh,self.xh,ad,ld)
            return np.block([[Lll,Llh],[Llh.T,Hhh]])+np.eye(len(self.y))*self.noise
        self.matrix=matrix
        def loss(theta):
            try:
                cf=cho_factor(matrix(theta),lower=True)
                return .5*self.y@cho_solve(cf,self.y)+np.log(np.diag(cf[0])).sum()+.5*len(self.y)*np.log(2*np.pi)
            except np.linalg.LinAlgError:return 1e30
        bounds=[(-3,3)]+[(-3,1)]*self.dim+[(-4,3)]+[(-3,1)]*self.dim+[(-3,3)]
        starts=[np.r_[0,np.full(self.dim,np.log(l)),np.log(.5),np.full(self.dim,np.log(l)),1.] for l in (.15,.4,.8)]
        attempts=[minimize(loss,s,method='L-BFGS-B',bounds=bounds,options={'maxiter':100,'ftol':1e-8}) for s in starts]
        valid=[]
        for attempt in attempts:
            if not np.isfinite(attempt.x).all() or not np.isfinite(attempt.fun) or attempt.fun>=1e29:continue
            try:
                K=matrix(attempt.x)
                if np.isfinite(K).all():cho_factor(K,lower=True);valid.append(attempt)
            except np.linalg.LinAlgError:pass
        if not valid:raise ArithmeticError('no finite SPD multi-fidelity GP optimizer iterate')
        converged=[a for a in valid if a.success];best=min(converged or valid,key=lambda a:a.fun);self.theta=best.x
        self.cf=cho_factor(matrix(best.x),lower=True);self.alpha=cho_solve(self.cf,self.y)
        self.fit_record={'theta':self.theta.tolist(),'rho':float(self.theta[-1]),'jitter_normalized_variance':self.noise,'covariance_condition_number':float(np.linalg.cond(matrix(best.x))),
            'mean_L':self.ml,'mean_H':self.mh,'output_scale':self.scale,'negative_log_likelihood':float(best.fun),
            'optimizer_attempts':[{'success':bool(a.success),'fun':float(a.fun) if np.isfinite(a.fun) else None,'message':str(a.message)} for a in attempts],
            'fit_status':'converged_MLE' if converged else 'finite_nonconverged_fallback',
            'fallback_policy':'prefer finite SPD converged starts; if none, lowest finite SPD iterate flagged, no test-based retries',
            'uncertainty_scope':'latent output conditional on fitted hyperparameters; no hyperparameter or natural-event uncertainty',
            'latent_models':'independent L and delta prior; full joint posterior cross-covariance retained'}
        return self

    def cross(self,x,level):
        al,ll,ad,ld,rho=self.unpack(self.theta)
        if level=='L':return np.c_[matern52(x,self.xl,al,ll),rho*matern52(x,self.xh,al,ll)]
        return np.c_[rho*matern52(x,self.xl,al,ll),rho*rho*matern52(x,self.xh,al,ll)+matern52(x,self.xh,ad,ld)]

    def prior(self,x,y,level_x='H',level_y='H'):
        al,ll,ad,ld,rho=self.unpack(self.theta)
        k=matern52(x,y,al,ll)
        if level_x==level_y=='L':return k
        if level_x!=level_y:return rho*k
        return rho*rho*k+matern52(x,y,ad,ld)

    def predict(self,x,level='H'):
        c=self.cross(x,level);mean=c@self.alpha
        var=np.diag(self.prior(x,x,level,level))-np.sum(c*cho_solve(self.cf,c.T).T,axis=1)
        if var.min() < -1e-7:raise ArithmeticError('negative GP posterior variance')
        return (self.ml if level=='L' else self.mh)+self.scale*mean,self.scale*np.sqrt(np.maximum(0,var))

    def posterior_cross(self,x,y,level_x='H',level_y='H'):
        return self.scale**2*(self.prior(x,y,level_x,level_y)-self.cross(x,level_x)@cho_solve(self.cf,self.cross(y,level_y).T))

    def acquisition(self,candidates,integration_grid,cost_low,cost_high,explore=False):
        m,s=self.predict(integration_grid)
        w=np.ones_like(m) if explore else np.exp(-.5*(m/np.maximum(s,1e-8))**2)
        w=w/max(w.sum(),1e-30)
        result={}
        for level,cost in [('L',cost_low),('H',cost_high)]:
            _,sq=self.predict(candidates,level)
            cov=self.posterior_cross(integration_grid,candidates,'H',level)
            result[level]=np.sum(w[:,None]*cov**2,axis=0)/(sq**2+self.noise*self.scale**2)/cost
        return result


def minmod3(a,b,c):
    return np.where((a*b>0)&(a*c>0),np.sign(a)*np.minimum(np.minimum(abs(a),abs(b)),abs(c)),0.)


def exact_wall_pressure(h,w):


    h=np.asarray(h,dtype=float);w=np.broadcast_to(np.asarray(w,dtype=float),h.shape)
    c=np.sqrt(G*h);cs=np.maximum(0,c+.5*w);hs=np.array(cs*cs/G,copy=True)
    positive=(w>0)&(h>0)
    if np.any(positive):
        f=w[positive]/c[positive];lo=np.ones_like(f);hi=1+2*f+f*f
        for _ in range(60):
            r=.5*(lo+hi);v=(r-1)*np.sqrt((r+1)/(2*r));lo=np.where(v<f,r,lo);hi=np.where(v>=f,r,hi)
        hs[positive]=h[positive]*.5*(lo+hi)
    return np.where(h>0,.5*G*hs*hs,0.)


def muscl_faces(V,Z,boundary):
    h=V[...,0];vel=np.divide(V[...,1:],h[...,None],out=np.zeros_like(V[...,1:]),where=h[...,None]>DRY)
    W=np.concatenate([(h+Z)[...,None],vel,Z[...,None]],axis=-1)
    sx=np.zeros_like(W);sy=np.zeros_like(W)
    sx[:,1:-1]=minmod3(2*(W[:,1:-1]-W[:,:-2]),.5*(W[:,2:]-W[:,:-2]),2*(W[:,2:]-W[:,1:-1]))
    sy[1:-1]=minmod3(2*(W[1:-1]-W[:-2]),.5*(W[2:]-W[:-2]),2*(W[2:]-W[1:-1]))

    badx=np.zeros_like(h,dtype=bool);bady=np.zeros_like(h,dtype=bool)
    badx[:,1:-1]=(h[:,1:-1]<=DRY)|(h[:,:-2]<=DRY)|(h[:,2:]<=DRY)
    bady[1:-1]=(h[1:-1]<=DRY)|(h[:-2]<=DRY)|(h[2:]<=DRY)
    badx|=h<.5*abs(sx[...,0]-sx[...,3]);bady|=h<.5*abs(sy[...,0]-sy[...,3])
    sx[badx]=0.;sy[bady]=0.
    if boundary=='slope':
        sx[:,0]=sx[:,1];sx[:,-1]=sx[:,-2];sy[0]=sy[1];sy[-1]=sy[-2]
    if boundary=='periodic':
        sx[:,0]=sx[:,-2];sx[:,-1]=sx[:,1];sy[0]=sy[-2];sy[-1]=sy[1]
    def state(A):
        hh=np.maximum(0,A[...,0]-A[...,3])
        return np.concatenate([hh[...,None],hh[...,None]*A[...,1:3]],axis=-1),A[...,3]
    left,zleft=state(W[1:-1,:-1]+.5*sx[1:-1,:-1]);right,zright=state(W[1:-1,1:]-.5*sx[1:-1,1:])
    top,ztop=state(W[:-1,1:-1]+.5*sy[:-1,1:-1]);bottom,zbottom=state(W[1:,1:-1]-.5*sy[1:,1:-1])
    wp=W[1:-1,1:-1];spx=sx[1:-1,1:-1];spy=sy[1:-1,1:-1]
    um,zm=state(wp-.5*spx);up,zp=state(wp+.5*spx)
    vm,ym=state(wp-.5*spy);vp,yp=state(wp+.5*spy)
    source_x=-.5*G*(um[...,0]+up[...,0])*(zp-zm)
    source_y=-.5*G*(vm[...,0]+vp[...,0])*(yp-ym)
    return left,right,zleft,zright,top,bottom,ztop,zbottom,source_x,source_y


def hll_pair(left, right, zl, zr, axis):

    hL, hR = left[..., 0], right[..., 0]
    uL = np.divide(left[..., 1:], hL[..., None], out=np.zeros_like(left[..., 1:]), where=hL[..., None] > DRY)
    uR = np.divide(right[..., 1:], hR[..., None], out=np.zeros_like(right[..., 1:]), where=hR[..., None] > DRY)
    zs = np.maximum(zl, zr)
    hl, hr = np.maximum(0, hL + zl - zs), np.maximum(0, hR + zr - zs)
    UL = np.concatenate([hl[..., None], hl[..., None]*uL], axis=-1)
    UR = np.concatenate([hr[..., None], hr[..., None]*uR], axis=-1)
    nl, nr = uL[..., axis], uR[..., axis]
    cl, cr = np.sqrt(G*hl), np.sqrt(G*hr)
    sl = np.minimum(nl-cl, nr-cr)
    sr = np.maximum(nl+cl, nr+cr)

    sl = np.where(hr <= 0., nl-cl, sl)
    sr = np.where(hr <= 0., nl+2*cl, sr)
    sl = np.where(hl <= 0., nr-2*cr, sl)
    sr = np.where(hl <= 0., nr+cr, sr)
    FL, FR = UL*nl[..., None], UR*nr[..., None]
    FL[..., axis+1] += .5*G*hl**2
    FR[..., axis+1] += .5*G*hr**2
    den = sr-sl
    middle = np.divide(sr[..., None]*FL-sl[..., None]*FR+(sl*sr)[..., None]*(UR-UL),
                       den[..., None], out=np.zeros_like(UL), where=den[..., None] > DRY)
    flux = np.where((sl >= 0)[..., None], FL, np.where((sr <= 0)[..., None], FR, middle))
    return flux, .5*G*(hL*hL-hl*hl), .5*G*(hR*hR-hr*hr)


def triangle_cdf(t, duration):
    r=np.clip(np.asarray(t)/duration,0,1)
    return np.where(r<=.5,2*r*r,1-2*(1-r)**2)


def pulse_volume(t0,t1,volume,duration,shape='triangle'):
    if shape=='double':

        half=duration/2
        return volume*.5*(triangle_cdf(t1,half)-triangle_cdf(t0,half)+
                          triangle_cdf(t1-half,half)-triangle_cdf(t0-half,half))
    return volume*(triangle_cdf(t1,duration)-triangle_cdf(t0,duration))


def euler_step(U,z,dx,dy,dt,mu=0.,xi=1e30,boundary='wall',inlet=None,t=0.):
    V=np.pad(U,((1,1),(1,1),(0,0)),mode='edge')
    Z=np.pad(z,((1,1),(1,1)),mode='edge')
    if boundary=='wall':
        V[:,0,1]*=-1;V[:,-1,1]*=-1;V[0,:,2]*=-1;V[-1,:,2]*=-1
    elif boundary=='slope':
        Z[:,0]=2*Z[:,1]-Z[:,2];Z[:,-1]=2*Z[:,-2]-Z[:,-3]
        Z[0]=2*Z[1]-Z[2];Z[-1]=2*Z[-2]-Z[-3]
    elif boundary=='periodic':
        V=np.pad(U,((1,1),(1,1),(0,0)),mode='wrap');Z=np.pad(z,((1,1),(1,1)),mode='wrap')
    elif boundary=='outflow':


        V[:,0]=0.;V[:,-1]=0.;V[0]=0.;V[-1]=0.
    elif boundary!='open':raise ValueError(boundary)
    if inlet is not None:


        V[1:-1,-1]=U[:,-1];V[1:-1,-1,1]*=-1
    if FAST is not None:
        fx,fy,xl,xr,yl,yr,source_x,source_y,left,zleft=FAST.fluxes(V,Z,ORDER,boundary,inlet is not None)
    else:
        source_x=np.zeros_like(z);source_y=np.zeros_like(z)
        if ORDER==2:
            left,right,zleft,zright,top,bottom,ztop,zbottom,source_x,source_y=muscl_faces(V,Z,boundary)
        else:
            left,right,zleft,zright=V[1:-1,:-1],V[1:-1,1:],Z[1:-1,:-1],Z[1:-1,1:]
            top,bottom,ztop,zbottom=V[:-1,1:-1],V[1:,1:-1],Z[:-1,1:-1],Z[1:,1:-1]


        if boundary=='wall':
            left[:,0]=right[:,0];left[:,0,1]*=-1;zleft[:,0]=zright[:,0]
            top[0]=bottom[0];top[0,:,2]*=-1;ztop[0]=zbottom[0]
            bottom[-1]=top[-1];bottom[-1,:,2]*=-1;zbottom[-1]=ztop[-1]
        if boundary=='wall' or inlet is not None:
            right[:,-1]=left[:,-1];right[:,-1,1]*=-1;zright[:,-1]=zleft[:,-1]
        fx,xl,xr=hll_pair(left,right,zleft,zright,0)
        fy,yl,yr=hll_pair(top,bottom,ztop,zbottom,1)
        if boundary=='wall':
            vx=np.divide(right[:,0,1],right[:,0,0],out=np.zeros(len(z)),where=right[:,0,0]>DRY)
            fx[:,0]=0.;fx[:,0,1]=exact_wall_pressure(right[:,0,0],-vx)
            vy=np.divide(bottom[0,:,2],bottom[0,:,0],out=np.zeros(z.shape[1]),where=bottom[0,:,0]>DRY)
            fy[0]=0.;fy[0,:,2]=exact_wall_pressure(bottom[0,:,0],-vy)
            vy=np.divide(top[-1,:,2],top[-1,:,0],out=np.zeros(z.shape[1]),where=top[-1,:,0]>DRY)
            fy[-1]=0.;fy[-1,:,2]=exact_wall_pressure(top[-1,:,0],vy)
        if boundary=='wall' or inlet is not None:
            vx=np.divide(left[:,-1,1],left[:,-1,0],out=np.zeros(len(z)),where=left[:,-1,0]>DRY)
            fx[:,-1]=0.;fx[:,-1,1]=exact_wall_pressure(left[:,-1,0],vx)
    vin=0.;fr=None;energy_in=0.;pressure_in=0.;inlet_old_flux=None
    if inlet is not None:
        weights=np.asarray(inlet.get('weights',inlet['rows']),dtype=float)
        mask=weights>0;speed=inlet['speed'];width=float(weights.sum()*dy)
        vin=float(pulse_volume(t,t+dt,inlet['volume'],inlet['duration'],inlet.get('shape','triangle')))
        q=vin/(dt*width);h=q/speed
        fr=speed/math.sqrt(G*h) if h>0 else None
        if h>0 and fr<=1:raise ArithmeticError('inlet not supercritical')

        inlet_old_flux=fx[mask,-1].copy()
        incoming=np.zeros_like(inlet_old_flux);incoming[:,0]=-q;incoming[:,1]=q*speed+.5*G*h*h
        hface2=[];wall_pressure=[]
        for j in np.flatnonzero(mask):
            bednodes=inlet.get('face_bed_samples',{}).get(int(j),np.array([zleft[j,-1]]))
            hh=np.maximum(0,left[j,-1,0]+zleft[j,-1]-bednodes)
            hface2.append(float(np.mean(hh**2)))
            if vin==0:
                vv=left[j,-1,1:]/left[j,-1,0] if left[j,-1,0]>DRY else np.zeros(2)
                wall_pressure.append(float(exact_wall_pressure(hh,vv[0]).mean()))
        hface2=np.array(hface2)
        if vin==0:


            incoming[:,1]=wall_pressure
        fx[mask,-1]=inlet_old_flux*(1-weights[mask,None])+incoming*weights[mask,None]
        pressure_in=dt*dy*float(np.sum((incoming[:,1]-q*speed)*weights[mask]))


        xl[mask,-1]=.5*G*(left[mask,-1,0]**2-hface2)*weights[mask]
        xr[mask,-1]=0.
        zin=inlet.get('bed_mean',float(np.average(z[mask,-1],weights=weights[mask])))
        energy_in=vin*(.5*speed**2+G*(zin-float(z.min())+h))
    bflux=dy*(fx[:,-1].sum(0)-fx[:,0].sum(0))+dx*(fy[-1].sum(0)-fy[0].sum(0))


    def head(A,B):
        h=A[...,0];v=np.divide(A[...,1:],h[...,None],out=np.zeros_like(A[...,1:]),where=h[...,None]>DRY)
        return .5*(v*v).sum(-1)+G*(B+h)
    H=head(U,z-float(z.min()))
    eb=dy*(np.sum(fx[:,-1,0]*H[:,-1])-np.sum(fx[:,0,0]*H[:,0]))+dx*(np.sum(fy[-1,:,0]*H[-1])-np.sum(fy[0,:,0]*H[0]))
    if inlet is not None:


        eb-=dy*float(np.sum(-q*weights[mask]*H[mask,-1]))
        eb-=energy_in/dt
    slope=np.zeros_like(U)
    slope[...,1]=-(xl[:,1:]-xr[:,:-1])/dx+source_x/dx
    slope[...,2]=-(yl[1:]-yr[:-1])/dy+source_y/dy
    result=U-dt/dx*(fx[:,1:]-fx[:,:-1])-dt/dy*(fy[1:]-fy[:-1])+dt*slope
    if result[...,0].min() < 0 or not np.isfinite(result).all():
        raise ArithmeticError('negative thickness or nonfinite state; no clipping')
    pre=result[...,1:].copy();h=result[...,0]
    vel=np.divide(pre,h[...,None],out=np.zeros_like(pre),where=h[...,None]>DRY)
    speed=np.linalg.norm(vel,axis=-1)
    r=np.maximum(0,speed-dt*mu*G)
    a=np.divide(dt*G,xi*h,out=np.zeros_like(h),where=h>DRY)
    newspeed=2*r/(1+np.sqrt(1+4*a*r))
    ratio=np.ones_like(speed)
    ratio[h>DRY]=np.divide(newspeed,speed,out=np.zeros_like(speed),where=speed>0)[h>DRY]
    result[...,1:]*=ratio[...,None]
    friction_impulse=(result[...,1:]-pre).sum((0,1))*dx*dy
    ke_before=.5*np.divide((pre*pre).sum(-1),h,out=np.zeros_like(h),where=h>0).sum()*dx*dy
    ke_after=.5*np.divide((result[...,1:]*result[...,1:]).sum(-1),h,out=np.zeros_like(h),where=h>0).sum()*dx*dy
    before_clear=result[...,1:].copy()
    result[...,1:][h<=DRY]=0
    clear=(result[...,1:]-before_clear).sum((0,1))*dx*dy
    ke_after_clear=.5*np.divide((result[...,1:]*result[...,1:]).sum(-1),h,out=np.zeros_like(h),where=h>0).sum()*dx*dy
    return result,{'boundary':dt*bflux,'bed':dt*slope[...,1:].sum((0,1))*dx*dy,
                   'friction':friction_impulse,'dry_clear':clear,'vin':vin,'Fr':fr,'energy_in':energy_in,
                   'energy_boundary':dt*eb,'friction_energy_loss':float(ke_before-ke_after),
                   'dry_clear_energy_loss':float(ke_after-ke_after_clear),
                   'boundary_faces':dt*np.concatenate([-dy*fx[:,0],dy*fx[:,-1],-dx*fy[0],dx*fy[-1]],axis=0),
                   'inlet_advective_impulse':vin*inlet['speed'] if inlet is not None else 0.,
                   'prescribed_inflow_pressure_impulse':pressure_in if vin>0 else 0.,
                   'closed_gate_pressure_impulse':pressure_in if vin==0 else 0.,
                   'inlet_pressure_impulse':pressure_in}


def step(U,z,dx,dy,dt,mu=0.,xi=1e30,boundary='wall',inlet=None,t=0.):
    one,b1=euler_step(U,z,dx,dy,dt,mu,xi,boundary,inlet,t)
    if ORDER==1:return one,b1
    two,b2=euler_step(one,z,dx,dy,dt,mu,xi,boundary,inlet,t)
    result=.5*U+.5*two
    budget={k:.5*(b1[k]+b2[k]) for k in b1 if k!='Fr'}


    h=result[...,0];dry=h<=DRY;before=result[...,1:].copy()
    extra_energy=.5*np.divide((before*before).sum(-1),h,out=np.zeros_like(h),where=(h>0)&dry).sum()*dx*dy
    result[...,1:][dry]=0.
    budget['dry_clear']+=(result[...,1:]-before).sum((0,1))*dx*dy
    budget['dry_clear_energy_loss']+=float(extra_energy)
    fs=[v for v in [b1['Fr'],b2['Fr']] if v is not None];budget['Fr']=min(fs) if fs else None
    return result,budget


def energy(U,z,area):
    h=U[...,0];v=np.divide(U[...,1:],h[...,None],out=np.zeros_like(U[...,1:]),where=h[...,None]>DRY)
    return float((.5*h*(v*v).sum(-1)+.5*G*h*h+G*h*z).sum()*area)


def evolve(U,z,dx,dy,end,mu=0.,xi=1e30,boundary='wall',cfl=.35,maxdt=1.,inlet=None,target=None,save_every=30.,trace_path=None,legal_exit_mask=None,exit_threshold_m3=None):
    started=time.perf_counter();t=0.;steps=0;area=dx*dy
    initial=U.sum((0,1))*area;boundary_sum=np.zeros(3);bed=np.zeros(2);friction=np.zeros(2);clear=np.zeros(2)
    zenergy=z-float(z.min())
    vin=0.;ein=0.;ebound=0.;edrag=0.;edry=0.;pinadv=0.;pinpressure=0.;einitial=energy(U,zenergy,area);frmin=math.inf;masserr=0.;momerr=0.;minh=float(U[...,0].min())
    ledgers=[];faces=[];vin_all=0.;vout_all=0.;max_energy=einitial;pin_flow=0.;pin_gate=0.
    side_out=np.zeros(4);illegal_out=0.;first_illegal=None
    if exit_threshold_m3 is None:exit_threshold_m3=1e-6*max(initial[0]+(inlet['volume'] if inlet is not None else 0),1)
    hmax=U[...,0].copy();vmax=np.zeros_like(z);event_vmax=np.zeros_like(z);snap=[];series=[];nextsave=0.
    while t<end-1e-10:
        h=U[...,0];vel=np.divide(U[...,1:],h[...,None],out=np.zeros_like(U[...,1:]),where=h[...,None]>DRY)
        c=np.sqrt(G*h);rate=float(np.max((abs(vel[...,0])+c)/dx+(abs(vel[...,1])+c)/dy))
        if inlet is not None:
            w=np.asarray(inlet.get('weights',inlet['rows']),dtype=float).sum()*dy
            hpeak=2*inlet['volume']/inlet['duration']/w/inlet['speed']
            cp=math.sqrt(G*hpeak)
            rate=max(rate,(inlet['speed']+cp)/dx+cp/dy)
        dt=min(maxdt,end-t,cfl/rate if rate>0 else maxdt)
        U,b=step(U,z,dx,dy,dt,mu,xi,boundary,inlet,t)
        t+=dt;steps+=1;boundary_sum+=b['boundary'];bed+=b['bed'];friction+=b['friction'];clear+=b['dry_clear']
        vin+=b['vin'];ein+=b['energy_in']
        ebound+=b['energy_boundary'];edrag+=b['friction_energy_loss']
        edry+=b['dry_clear_energy_loss']
        pinadv+=b['inlet_advective_impulse'];pinpressure+=b['inlet_pressure_impulse']
        pin_flow+=b['prescribed_inflow_pressure_impulse'];pin_gate+=b['closed_gate_pressure_impulse']


        vin_all+=float(-np.minimum(b['boundary_faces'][:,0],0).sum());vout_all+=float(np.maximum(b['boundary_faces'][:,0],0).sum())
        if b['Fr'] is not None:frmin=min(frmin,b['Fr'])
        total=U.sum((0,1))*area
        bf=b['boundary_faces'];cuts=[0,len(z),2*len(z),2*len(z)+z.shape[1],2*len(z)+2*z.shape[1]]
        for k in range(4):side_out[k]+=float(np.maximum(bf[cuts[k]:cuts[k+1],0],0).sum())
        if legal_exit_mask is not None:
            illegal_out+=float(np.maximum(bf[~legal_exit_mask,0],0).sum())
            if first_illegal is None and illegal_out>exit_threshold_m3:first_illegal=t
        et=energy(U,zenergy,area);max_energy=max(max_energy,et)
        if trace_path is not None:
            faces.append(b['boundary_faces'])
            ledgers.append([t,dt,*total,*b['bed'],*b['friction'],*b['dry_clear'],et,b['energy_boundary'],b['friction_energy_loss'],b['dry_clear_energy_loss'],b['vin'],b['Fr'] or 0.,b['inlet_advective_impulse'],b['prescribed_inflow_pressure_impulse'],b['closed_gate_pressure_impulse'],b['energy_in']])
        masserr=max(masserr,abs(total[0]-initial[0]+boundary_sum[0]))
        momerr=max(momerr,float(np.max(abs(total[1:]-initial[1:]+boundary_sum[1:]-bed-friction-clear))))
        minh=min(minh,float(U[...,0].min()));hmax=np.maximum(hmax,U[...,0])
        v=np.divide(U[...,1:],U[...,0,None],out=np.zeros_like(U[...,1:]),where=U[...,0,None]>DRY)
        vmax=np.maximum(vmax,np.linalg.norm(v,axis=-1))
        event_vmax=np.maximum(event_vmax,np.where(U[...,0]>=.1,np.linalg.norm(v,axis=-1),0.))
        target_q=0.;peak_east=0.;peak_north=0.;peak_tangent=0.
        if isinstance(target,dict):
            line=sum(w*(U[j,:,:] if target.get('axis')=='y' else U[:,j,:]) for w,j in zip(target['weights'],target['columns']))
            vals=line[:,0]
            target_h=float(vals[target['rows']].max())
            active=np.flatnonzero(target['rows']);k=int(active[np.argmax(vals[active])])
            if vals[k]>DRY:peak_east=float(line[k,1]/vals[k]);peak_north=float(-line[k,2]/vals[k])
            lengths=target.get('quadrature_lengths',target['rows']*(dx if target.get('axis')=='y' else dy))
            normal=2 if target.get('axis')=='y' else 1
            target_q=float(np.sum(line[:,normal]*lengths)*target.get('normal_sign',1))
            peak_tangent=float(target.get('tangent_coordinates',np.arange(len(line)))[k])
        else:target_h=float(U[...,0][target].max()) if target is not None else 0.
        series.append((t,target_h,float(total[0]),float(np.linalg.norm(U[...,1:],axis=-1).sum()*area),target_q,peak_east,peak_north,peak_tangent))
        if t>=nextsave or t>=end-1e-10:
            snap.append((t,U.copy()));nextsave+=save_every
        if steps>100000:raise RuntimeError('step budget exceeded')
    total=U.sum((0,1))*area
    diag={'steps':steps,'wall_s':time.perf_counter()-started,'min_h':minh,'mass_error_abs_m3':masserr,
          'mass_error_rel':masserr/max(initial[0]+vin,1),'momentum_error_abs':momerr,
          'momentum_error_rel':momerr/max(float(np.linalg.norm(bed))+float(np.linalg.norm(boundary_sum[1:])),1),
          'initial':initial.tolist(),'final':total.tolist(),'boundary_integral':boundary_sum.tolist(),
          'bed_impulse':bed.tolist(),'friction_impulse':friction.tolist(),'dry_clear_impulse':clear.tolist(),
          'inlet_volume_m3':vin,'min_inlet_Fr':frmin if math.isfinite(frmin) else None,
          'total_inflow_m3':vin_all,'total_outflow_m3':vout_all,'dry_threshold_m':DRY,
          'outflow_by_side_m3':dict(zip(['west','east','north','south'],side_out.tolist())),
          'illegal_outflow_m3':illegal_out,'first_significant_illegal_outflow_s':first_illegal,'illegal_outflow_threshold_m3':exit_threshold_m3,
          'whole_window_valid':first_illegal is None,'valid_until_s':first_illegal if first_illegal is not None else end,
          'inlet_advective_impulse_west_per_density':pinadv,'inlet_pressure_impulse_west_per_density':pinpressure,
          'prescribed_inflow_pressure_impulse_west_per_density':pin_flow,'closed_gate_pressure_impulse_west_per_density':pin_gate,
          'legacy_inlet_pressure_definition':'total inlet-strip pressure including post-pulse closed gate; not solely prescribed inflow',
          'energy_in_per_density':ein,'energy_final_per_density':energy(U,zenergy,area),'energy_reference_elevation_m':float(z.min()),
          'energy_initial_per_density':einitial,'energy_boundary_outward_per_density':ebound,
          'friction_energy_loss_per_density':edrag,
          'dry_clear_energy_loss_per_density':edry,'maximum_energy_per_density':max_energy,
          'energy_residual_per_density':energy(U,zenergy,area)-einitial+ebound+edrag+edry,
          'energy_status':'first-order boundary head flux; residual combines transport dissipation, bed discretization and source splitting; no exact entropy claim',
          'independent_review':'pending'}
    if trace_path is not None:
        np.savez_compressed(trace_path,ledger=ledgers,boundary_faces=faces,initial=initial,initial_energy=einitial,
            columns='t,dt,V,Px,Psouth,bed_x,bed_south,drag_x,drag_south,clear_x,clear_south,E,Eboundary,Eloss_drag,Eloss_clear,Vin,Fr,Pin_adv,Pin_flow_pressure,Pin_gate_pressure,Ein',
            face_order='west,east,north,south; each outward and time-integrated; momentum per density',dx=dx,dy=dy)
    return U,diag,np.array(series),snap,hmax,vmax,event_vmax


def fp(p):
    import hashlib
    p=Path(p);return {'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}


def write_json(p,data):
    Path(p).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def runtime_versions():
    from importlib.metadata import version,PackageNotFoundError
    result={}
    for name in ('numpy','scipy','matplotlib','pillow','rasterio','pyproj'):
        try:result[name]=version(name)
        except PackageNotFoundError:result[name]=None
    return result


def benchmarks(out):
    checks=[]
    def add(name,ok,**data):
        row={'test':name,'passed':bool(ok),**data};checks.append(row);print(json.dumps(row),flush=True)
    time_grid=np.linspace(0,960,20001);dt=time_grid[1]-time_grid[0]
    single=np.asarray(pulse_volume(time_grid[:-1],time_grid[1:],8e5,960,'triangle'))/dt
    double=np.asarray(pulse_volume(time_grid[:-1],time_grid[1:],8e5,960,'double'))/dt
    volume_error=float(abs(single.sum()-double.sum())*dt/8e5)
    pressure_error=float(abs(np.sum(single**2)-np.sum(double**2))/np.sum(single**2))
    add('equal_volume_support_peak_double_pulse',volume_error<1e-10 and pressure_error<1e-6 and abs(single.max()-double.max())/single.max()<1e-3,
        volume_relative_difference=volume_error,squared_flux_integral_relative_difference=pressure_error,total_support_s=960,design='two V/2,T/2 triangles adjacent in [0,T]')
    x,y=np.meshgrid(np.linspace(-2,2,41),np.linspace(-2,2,37));z=.4*np.cos(x)*np.cos(y)+.08*x+.1*y
    U=np.zeros((*z.shape,3));U[...,0]=np.maximum(0,.25-z);old=U.copy()
    U,d,*_=evolve(U,z,1,1,10,maxdt=.1)
    dh=float(abs(U[...,0]-old[...,0]).max());speed=float(np.linalg.norm(np.divide(U[...,1:],U[...,0,None],out=np.zeros_like(U[...,1:]),where=U[...,0,None]>DRY),axis=-1).max())
    add('lake_wet_dry',dh<1e-9 and speed<1e-8,depth_error=dh,speed=speed,diagnostics=d)
    errs=[]
    for dx in (10.,5.,2.5):
        x=np.arange(-1000+dx/2,1000,dx);z=np.zeros((1,len(x)));U=np.zeros((*z.shape,3));U[0,:,0]=np.where(x<0,10,0)
        U,d,*_=evolve(U,z,dx,1e8,20,boundary='open')
        c=math.sqrt(G*10);r=x/20;truth=np.where(r<-c,10,np.where(r<=2*c,(2*c-r)**2/(9*G),0))
        err=float(abs(U[0,:,0]-truth).sum()/truth.sum());errs.append(err)
        np.savez(out/f'ritter_{dx:g}.npz',x=x,h=U[0,:,0],exact=truth)
        add(f'ritter_{dx:g}',err<.15,error=err,diagnostics=d)
    add('ritter_convergence',all(a>b for a,b in zip(errs,errs[1:])) and errs[-1]<.08,errors=errs)
    xx,yy=np.meshgrid(np.arange(-20,21),np.arange(-20,21));z=np.zeros(xx.shape);U=np.zeros((*z.shape,3));U[...,0]=np.where(xx*xx+yy*yy<25,2,0)
    U,d,*_=evolve(U,z,1,1,3,trace_path=out/'dry_spread_budget.npz')
    closed_mass_error=abs(d['final'][0]-d['initial'][0])/d['initial'][0]
    add('dry_spread',d['mass_error_rel']<1e-9 and d['min_h']>=0 and closed_mass_error<1e-10 and d['total_inflow_m3']<1e-10 and d['total_outflow_m3']<1e-10,
        diagnostics=d,closed_storage_relative_error=closed_mass_error,symmetry=float(abs(U[...,0]-U[...,0].T).max()))
    add('closed_flat_no_energy_creation',d['maximum_energy_per_density']<=d['energy_initial_per_density']*(1+1e-8),initial_energy=d['energy_initial_per_density'],maximum_energy=d['maximum_energy_per_density'])

    slope=.08;mu=.03;xi=500.;depth=2.;T=4.;A=G*(slope-mu);B=G/(xi*depth)
    exact=math.sqrt(A/B)*math.tanh(math.sqrt(A*B)*T);errs=[]
    for dx in (.5,.25,.125):
        z=np.broadcast_to(-slope*np.arange(16)*dx,(5,16)).copy();U=np.zeros((*z.shape,3));U[...,0]=depth
        U,d,*_=evolve(U,z,dx,dx,T,mu,xi,'slope',maxdt=dx*.02)
        v=float((U[...,1]/U[...,0]).mean());err=abs(v-exact)/exact;errs.append(err)
        add(f'slope_integrated_{dx:g}',err<.1,speed=v,exact=exact,error=err,diagnostics=d)
    add('slope_convergence',all(a>b for a,b in zip(errs,errs[1:])) and errs[-1]<.03,errors=errs)
    from scipy.integrate import solve_ivp
    for label,v0 in [('terminal',math.sqrt(A/B)),('superterminal',15.)]:
        dx=.125;z=np.broadcast_to(-slope*np.arange(16)*dx,(5,16)).copy();U=np.zeros((*z.shape,3));U[...,0]=depth;U[...,1]=depth*v0
        U,d,*_=evolve(U,z,dx,dx,T,mu,xi,'slope',maxdt=.0025)
        reference=float(solve_ivp(lambda t,v:A-B*v*v,[0,T],[v0],rtol=1e-11,atol=1e-12,method='DOP853').y[0,-1])
        observed=float((U[...,1]/U[...,0]).mean());err=abs(observed-reference)/reference
        add('slope_'+label,err<.01,reference=reference,observed=observed,error=err,diagnostics=d)
    z=np.broadcast_to(-.01*np.arange(20),(5,20)).copy();U=np.zeros((*z.shape,3));U[...,0]=2
    U,d,*_=evolve(U,z,1,1,2,.1,500,'slope',maxdt=.02)
    add('subthreshold_no_start',float(abs(U[...,1:]).max())<1e-12,max_momentum=float(abs(U[...,1:]).max()))
    z=np.zeros((8,40));U=np.zeros((*z.shape,3));U[...,0]=2;U[...,1]=.02
    U,d,*_=evolve(U,z,1,1,2,.1,500,'open',maxdt=.02)
    add('stop_no_reverse',float(U[...,1].min())>=0 and float(U[...,1].max())<1e-12,min_momentum=float(U[...,1].min()),max_momentum=float(U[...,1].max()))
    z=np.zeros((10,80));U=np.zeros((*z.shape,3));mask=np.ones(10,dtype=bool)
    inlet={'rows':mask,'speed':15.,'volume':1000.,'duration':20.}
    U,d,*_=evolve(U,z,5,5,60,.01,500,'outflow',inlet=inlet,trace_path=out/'inlet_budget.npz')
    add('inlet_open_budget',abs(d['inlet_volume_m3']-1000)<1e-3 and abs(d['total_inflow_m3']-1000)<1e-3 and d['mass_error_rel']<1e-8 and d['momentum_error_rel']<1e-8,diagnostics=d)

    z=np.zeros((12,40));U=np.zeros((*z.shape,3));U[:]=[2.,4.,1.5]
    old=U.copy();U,d,*_=evolve(U,z,3,3,2,boundary='open',maxdt=.05)
    add('oblique_open_uniform',float(abs(U-old).max())<1e-10 and d['momentum_error_rel']<1e-8,diagnostics=d)

    z=np.zeros((12,30));U=np.zeros((*z.shape,3));U[3:9,18:26,0]=1.;U[3:9,18:26,1]=15.;U[3:9,18:26,2]=2.
    U,d,*_=evolve(U,z,2,2,10,boundary='outflow',maxdt=.01,trace_path=out/'dry_exterior_budget.npz')
    add('dry_exterior_actual_exit',d['total_outflow_m3']>.8*d['initial'][0] and d['total_inflow_m3']<1e-9 and d['mass_error_rel']<1e-8 and d['momentum_error_rel']<1e-8,diagnostics=d,
        analytic_bound='initial 192m3, advective edge clears east boundary in <1.6s at15m/s; gravity spreading allowed, require >80% exits in10s')

    z=np.zeros((1,100));U=np.zeros((*z.shape,3));U[...,0]=1.
    U,d,*_=evolve(U,z,1,1e8,.5,boundary='outflow',maxdt=.01,trace_path=out/'subcritical_drain_budget.npz')
    add('subcritical_dry_exterior_drain',d['total_outflow_m3']>0 and d['total_inflow_m3']<1e-8 and d['mass_error_rel']<1e-8,diagnostics=d,
        restriction='wet interior at rest, dry exterior dam-break; no imposed downstream stage or backwater claim')
    U=np.zeros((5,10,3));U[...,0]=1.;z=np.zeros((5,10));A,_=step(U,z,1,1,.01,boundary='outflow');U[...,1]=1e-12;B,_=step(U,z,1,1,.01,boundary='outflow')
    jump=float(abs(A-B).max());add('boundary_small_perturbation_continuity',jump<1e-8,max_state_difference=jump,perturbation=1e-12)
    U=np.zeros((4,4,3));U[...,0]=DRY*.5;U[...,1]=1e-4;z=np.zeros((4,4));A,b=step(U,z,1,1,.01,boundary='open')
    add('dry_clear_separate_from_drag',float(abs(b['friction']).max())<1e-15 and abs(b['dry_clear'][0]+.0016)<1e-14,
        friction_impulse=b['friction'].tolist(),dry_clear_impulse=b['dry_clear'].tolist(),dry_clear_energy=b['dry_clear_energy_loss'])
    errs=[];T=.5;mu=.02;xi=300.;depth=2.;v0=20.;aa=G*mu;bb=G/(xi*depth)
    ref=math.sqrt(aa/bb)*math.tan(math.atan(v0*math.sqrt(bb/aa))-T*math.sqrt(aa*bb))
    for dt in (.05,.025,.0125):
        z=np.zeros((5,10));U=np.zeros((*z.shape,3));U[...,0]=depth;U[...,1]=depth*v0*.6;U[...,2]=depth*v0*.8
        U,d,*_=evolve(U,z,100,100,T,mu,xi,'open',maxdt=dt)
        obs=float(np.linalg.norm(U[...,1:]/U[...,0,None],axis=-1).mean());errs.append(abs(obs-ref)/ref)
    add('oblique_pure_drag_convergence',errs[-1]<.005 and all(a>b for a,b in zip(errs,errs[1:])),errors=errs,reference=ref)
    if ORDER==2:
        errors=[]
        for dx in (60.,30.,15.):
            ss=17/30;depth=1.85;dt=1e-4;z=np.broadcast_to(-ss*np.arange(16)*dx,(5,16)).copy();U=np.zeros((*z.shape,3));U[...,0]=depth
            A,b=step(U,z,dx,dx,dt,.05,500,'slope');observed=float((A[...,1]/A[...,0]).mean())/dt;reference=G*(ss-.05)
            err=abs(observed-reference)/reference;errors.append(err)
            add('thin_steep_start_'+str(dx),err<1e-4,acceleration=observed,reference=reference,relative_error=err)
        U=np.zeros((3,8,3));U[...,0]=1.;U[...,1]=np.array([12,12,12,12,12,12,11,10])[None,:];z=np.zeros((3,8))
        A,b=step(U,z,1,1,.001,boundary='wall')
        add('reconstructed_wall_gradient_no_water',abs(float(A[...,0].sum()-U[...,0].sum()))<1e-12 and float(abs(b['boundary_faces'][:,0]).max())<1e-12,
            volume_change=float(A[...,0].sum()-U[...,0].sum()),maximum_wall_mass_flux_integral=float(abs(b['boundary_faces'][:,0]).max()),wall_energy=b['energy_boundary'])
        weights=np.array([.2,.3,.5]);inlet={'rows':weights>0,'weights':weights,'speed':15.,'volume':100.,'duration':20.}
        A,b=step(U,z,1,1,.001,boundary='wall',inlet=inlet,t=30.)
        add('partial_gate_gradient_no_water',abs(float(A[...,0].sum()-U[...,0].sum()))<1e-12 and float(abs(b['boundary_faces'][:,0]).max())<1e-12,
            volume_change=float(A[...,0].sum()-U[...,0].sum()),maximum_wall_mass_flux_integral=float(abs(b['boundary_faces'][:,0]).max()))
        from scipy.optimize import brentq
        pressures=[];references=[]
        for w in (-15.,-math.sqrt(G),0.,2.,15.):
            p=float(exact_wall_pressure(1.,w));pressures.append(p)
            if w<=0:ref=.5/G*max(0,math.sqrt(G)+w/2)**4
            else:
                rr=brentq(lambda r:(r-1)*math.sqrt((r+1)/(2*r))-w/math.sqrt(G),1,1+2*w/math.sqrt(G)+(w/math.sqrt(G))**2,xtol=1e-13)
                ref=.5*G*rr*rr
            references.append(ref)
        add('exact_wall_rarefaction_vacuum_shock',max(abs(a-b) for a,b in zip(pressures,references))<1e-9 and pressures[0]==0,
            pressure=pressures,reference=references)
        dry_pressures=exact_wall_pressure(np.zeros(3),np.array([-15.,0.,15.]))
        add('exact_wall_dry_zero_pressure',float(abs(dry_pressures).max())==0,pressure=dry_pressures.tolist())
        U=np.zeros((3,8,3));U[...,0]=1.;U[...,1]=15.;z=np.zeros((3,8))
        inlet={'rows':weights>0,'weights':weights,'speed':15.,'volume':100.,'duration':20.,'face_bed_samples':{j:np.array([0.,2.]) for j in range(3)}}
        _,bb=euler_step(U,z,1,1,1e-6,boundary='wall',inlet=inlet,t=30.)
        expected=.5*float(exact_wall_pressure(1.,15.))*weights.sum()*1e-6
        add('partial_closed_dry_face_nodes',abs(bb['inlet_pressure_impulse']-expected)<1e-12,pressure_impulse=bb['inlet_pressure_impulse'],reference=expected)
    write_json(out/'CHECKS.json',{'checks':checks,'all_passed':all(c['passed'] for c in checks),'independent_review':'pending'})


def smooth_benchmark(out):
    errors=[]
    for n in (100,200,400):
        length=1000.;dx=length/n;x=(np.arange(n)+.5)*dx;c0=math.sqrt(2*G);amp=.1*c0;k=2*np.pi/length;T=10.
        c=c0+amp*np.sin(k*x);h=c*c/G;u=2*(c-c0)
        U=np.zeros((1,n,3));U[0,:,0]=h;U[0,:,1]=h*u;z=np.zeros((1,n))
        U,d,*_=evolve(U,z,dx,1e8,T,boundary='periodic',maxdt=10.)
        exact_c=np.full_like(x,c0)
        for _ in range(100):exact_c=c0+amp*np.sin(k*(x-(3*exact_c-2*c0)*T))
        exact_h=exact_c**2/G;err=float(abs(U[0,:,0]-exact_h).mean());errors.append(err)
        np.savez(out/f'smooth_{n}.npz',x=x,h=U[0,:,0],exact_h=exact_h)
        print(json.dumps({'n':n,'L1_h':err,'steps':d['steps']}),flush=True)
    orders=np.log2(np.array(errors[:-1])/np.array(errors[1:]))
    write_json(out/'SMOOTH_CHECK.json',{'errors':errors,'observed_orders':orders.tolist(),'passed':bool(all(a>b for a,b in zip(errors,errors[1:])) and orders[-1]>1.5),'scope':'wet smooth friction-free simple wave only; full friction integration is mixed order','reconstruction_order':ORDER})


def terrain(out,root,full=False):
    import csv,rasterio
    from rasterio.warp import reproject,Resampling
    from rasterio.transform import from_origin
    from pyproj import Transformer
    srcdir=root/'补充数据/吉隆数据 2/09_补充DEM'
    sources={'copernicus':srcdir/'Copernicus_DSM_COG_10_N28_00_E085_00_DEM.tif','aster':srcdir/'ASTGTMV003_N28E085/ASTGTMV003_N28E085_dem.tif'}

    bounds=[339840,3128160,350640,3136080] if full else [348000,3134520,350640,3135840]
    res=15.;transform=from_origin(bounds[0],bounds[3],res,res)
    records={}
    for name,p in sources.items():
        a=np.full((int((bounds[3]-bounds[1])/res),int((bounds[2]-bounds[0])/res)),np.nan)
        with rasterio.open(p) as s:
            reproject(rasterio.band(s,1),a,src_transform=s.transform,src_crs=s.crs,dst_transform=transform,dst_crs='EPSG:32645',resampling=Resampling.bilinear,dst_nodata=np.nan)
            records[name]={'source':fp(p),'original_crs':str(s.crs),'original_transform':list(s.transform),'nodata':s.nodata,'vertical_datum':'product documentation pending; not from geographic CRS','missing':int((~np.isfinite(a)).sum())}
        if not np.isfinite(a).all():raise ValueError('terrain missing values')
        np.savez_compressed(out/f'{name}_support.npz',z=a,bounds=bounds,resolution=res,crs='EPSG:32645')
    rows=list(csv.DictReader((root/'outputs/Q3_implementation_20261030/data_diagnostics/gyirong_author_mapped_path.csv').open()))
    tr=Transformer.from_crs(4326,32645,always_xy=True)
    coords=np.array([tr.transform(float(r['longitude']),float(r['latitude'])) for r in rows])
    np.savetxt(out/'mapped_path_utm.csv',coords[239:721] if full else coords[239:401],delimiter=',',header='easting_m,northing_m',comments='')
    write_json(out/'TERRAIN.json',{'sources':records,'bounds':bounds,'support_resolution_m':15,'flow_grid_candidates_m':[60,30],
        'entry':{'x_m':350640,'y_interval_m':[3135025,3135055],'width_m':30,'speed_m_s':30,'direction':'west','identity':'D_ONLY boundary near author node270, ~928m after C; no measured state'},
        'target':{'y_m':3129103.736725446,'x_range_m':[340500,341100],'identity':'port-reference cross-section near author699; no observed hydrograph'} if full else {'x_m':348720,'y_range_m':[3134880,3135600],'identity':'nested cross-section; not port'},
        'no_filling_or_elevation_correction':True,'hydraulic_qualification':'pending','full_event_reconstruction':False})


def sample_terrain(p,dx):
    from scipy.ndimage import map_coordinates
    f=np.load(p);a=f['z'];bounds=f['bounds'];res=float(f['resolution'])
    if str(f['crs'])!='EPSG:32645' or res<=0 or not np.isfinite(a).all() or a.ndim!=2:raise ValueError('finite EPSG32645 metre terrain support required')
    expected=((bounds[3]-bounds[1])/res,(bounds[2]-bounds[0])/res)
    if not np.allclose(expected,a.shape,rtol=0,atol=1e-9):raise ValueError('terrain bounds/resolution/shape mismatch')
    if not np.isfinite(dx) or dx<=0 or not np.allclose(np.array(expected)*res/dx,np.round(np.array(expected)*res/dx),rtol=0,atol=1e-9):raise ValueError('flow resolution must tile the fixed physical domain exactly')
    nx=int(round((bounds[2]-bounds[0])/dx));ny=int(round((bounds[3]-bounds[1])/dx))

    xx=bounds[0]+(np.arange(nx)+.5)*dx;yy=bounds[3]-(np.arange(ny)+.5)*dx
    cc,rr=np.meshgrid((xx-bounds[0])/res-.5,(bounds[3]-yy)/res-.5)
    return map_coordinates(a,[rr,cc],order=1,mode='nearest'),xx,yy


def pilot(out,inputs,dx,args):
    z,x,y=sample_terrain(inputs/f'{args.terrain_name}_support.npz',dx);U=np.zeros((*z.shape,3))
    width=float(args.entry_width);entry_lo=3135040-width/2;entry_hi=3135040+width/2
    if not np.isfinite(width) or width<=0 or entry_lo<y.min()-dx/2 or entry_hi>y.max()+dx/2:raise ValueError('entry width outside domain')
    weights=np.maximum(0,np.minimum(y+dx/2,entry_hi)-np.maximum(y-dx/2,entry_lo))/dx
    mask=weights>0
    from scipy.ndimage import map_coordinates
    f=np.load(inputs/f'{args.terrain_name}_support.npz');ys=np.linspace(entry_lo,entry_hi,121);bounds=f['bounds']
    zs=map_coordinates(f['z'],[(bounds[3]-ys)/15-.5,np.full(len(ys),f['z'].shape[1]-.5)],order=1,mode='nearest')
    zin=float(np.trapezoid(zs,ys)/width)
    inlet={'rows':mask,'weights':weights,'speed':args.speed,'volume':args.volume,'duration':args.duration,'shape':args.shape,'bed_mean':zin}
    inlet['face_bed_samples']={int(j):np.interp(np.linspace(max(y[j]-dx/2,entry_lo),min(y[j]+dx/2,entry_hi),61),ys,zs) for j in np.flatnonzero(mask)}
    write_json(out/'ENTRY_FACE.json',{'x_m':350640,'y_interval_m':[entry_lo,entry_hi],'continuous_bed_y_m':ys.tolist(),'continuous_bed_z_m':zs.tolist(),'bed_mean_m':zin,
        'grid_faces':[{'row':int(j),'center_y_m':float(y[j]),'fraction':float(weights[j]),'mapped_cell_bed_m':float(z[j,-1]),'face_mean_bed_m':float(inlet['face_bed_samples'][int(j)].mean()),'bed_mapping_difference_m':float(z[j,-1]-inlet['face_bed_samples'][int(j)].mean())} for j in np.flatnonzero(mask)],
        'normal_in_utm':[-1,0],'width_m':width,'speed_m_s':args.speed,'peak_physical_h_m':2*args.volume/args.duration/width/args.speed,'subcell_pressure':'fraction*g*h_physical^2/2, complementary face wall; half-cell bed pressure correction',
        'boundary_restriction':'dry exterior except upstream east gate; inlet closed outside pulse support','vector_convention':'internal [east,south]; exported north velocity=-internal south'})

    if args.full:
        yt=3129103.736725446;row=int(np.searchsorted(-y,-yt));alpha=(yt-y[row-1])/(y[row]-y[row-1])
        target={'axis':'y','columns':[row-1,row],'weights':[1-alpha,alpha],'rows':(x>=340500)&(x<=341100)}
        target.update({'quadrature_lengths':np.maximum(0,np.minimum(x+dx/2,341100)-np.maximum(x-dx/2,340500)),'normal_sign':1,'tangent_coordinates':x})
    else:
        col=int(np.searchsorted(x,348720));alpha=(348720-x[col-1])/(x[col]-x[col-1])
        target={'axis':'x','columns':[col-1,col],'weights':[1-alpha,alpha],'rows':(y>=3134880)&(y<=3135600)}
        target.update({'quadrature_lengths':np.maximum(0,np.minimum(y+dx/2,3135600)-np.maximum(y-dx/2,3134880)),'normal_sign':-1,'tangent_coordinates':y})
    legal=np.zeros(2*len(y)+2*len(x),dtype=bool)
    if args.full:legal[2*len(y)+len(x):]=(x>=340500)&(x<=341100)
    else:legal[:len(y)]=(y>=3134880)&(y<=3135600)
    U,d,series,snaps,hm,vm,event_vm=evolve(U,z,dx,dx,args.horizon,args.mu,args.xi,'outflow',maxdt=args.maxdt,cfl=args.cfl,inlet=inlet,target=target,trace_path=out/f'pilot_{dx:g}_budget.npz',legal_exit_mask=legal)
    np.savez_compressed(out/f'pilot_{dx:g}.npz',z=z,x=x,y=y,final=U,max_h=hm,max_speed=vm,max_speed_when_h_ge_0p1=event_vm,series=series,
        series_columns='time_s,target_max_h_m,volume_m3,total_momentum_norm_per_density,Q_downstream_m3_s,u_east_at_hpeak_m_s,u_north_at_hpeak_m_s,tangent_coordinate_at_hpeak_m',snapshot_t=[s[0] for s in snaps],snapshot_U=[s[1] for s in snaps])
    arrivals={}
    for threshold in (.05,.1,.2):
        hits=np.flatnonzero(series[:,1]>=threshold)
        prefix_valid=bool(len(hits) and series[hits[0],0]<=d['valid_until_s'])
        arrivals[str(threshold)]={'status':'arrived' if prefix_valid else ('valid_no_arrival_by_horizon' if d['whole_window_valid'] else 'unknown'),
            'time_bracket_s':[float(series[hits[0]-1,0]) if hits[0]>0 else 0.,float(series[hits[0],0])] if prefix_valid else None}
    finalspeed=np.linalg.norm(np.divide(U[...,1:],U[...,0,None],out=np.zeros_like(U[...,1:]),where=U[...,0,None]>DRY),axis=-1)
    d.update({'dx_m':dx,'shape':list(z.shape),'target_max_h':float(series[:,1].max()),'threshold_arrivals':arrivals,'status':arrivals['0.1']['status'],
        'run_termination':'out_of_domain' if not d['whole_window_valid'] else 'time_horizon','final_motion_state':'stopped' if finalspeed.max()<1e-3 else 'moving',
        'terrain_qualification':'pending','full_reach':args.full,'legal_exit_geometry':{'side':'south' if args.full else 'west','interval_m':[340500,341100] if args.full else [3134880,3135600]},
        'arrival_support':{'axis':target['axis'],'coordinate_m':3129103.736725446 if args.full else 348720,'interpolation_columns':target['columns'],'interpolation_weights':target['weights']},
        'field_convention':'saved U=(h,hu_east,hu_south); north momentum=-hu_south; lengths m, times s, velocity m/s'})
    d.update({'peak_h_anywhere_m':float(hm.max()),'maximum_speed_m_s':float(vm.max()),'inundated_area_m2':float((hm>=.1).sum()*dx*dx),
              'pulse_definition':'two adjacent half-volume/half-duration triangles within [0,T]' if args.shape=='double' else 'single triangle in [0,T]'})
    peak=int(np.argmax(series[:,1]));d.update({'target_Q_downstream_max_m3_s':float(series[:,4].max()),'target_Q_downstream_min_m3_s':float(series[:,4].min()),
        'target_peak_h_paired':dict(zip(['time_s','h_m','u_east_m_s','u_north_m_s','tangent_coordinate_m'],[float(series[peak,j]) for j in [0,1,5,6,7]])),
        'speed_envelope_for_figures':'max speed computed only at steps with local event h>=0.1m; raw max_speed separately retained','loads_not_computed':True})
    if series[peak,1]<.1:
        d['target_peak_h_paired']={'time_s':None,'h_m':float(series[peak,1]),'u_east_m_s':None,'u_north_m_s':None,'tangent_coordinate_m':None}
    d['ever_threshold_coverage_area_m2']=d['inundated_area_m2'];d['final_threshold_inundated_area_m2']=float((U[...,0]>=.1).sum()*dx*dx)
    d['arrival_support'].update({'crs':'EPSG:32645','tangent_interval_m':[340500,341100] if args.full else [3134880,3135600],'criterion':'instantaneous maximum vertical event depth >= threshold on fixed interpolated line'})
    write_json(out/f'pilot_{dx:g}.json',d);print(json.dumps(d),flush=True)


def run_morris(out,args):

    import csv,subprocess,sys,os,hashlib
    from concurrent.futures import ThreadPoolExecutor,as_completed
    names=['volume_m3','duration_s','mu','xi_m_s2']
    limits=np.array([[2e5,8e5],[720.,1200.],[.03,.15],[250.,1000.]])
    points,orders=morris_design(4,8,6,83031);points=np.round(points,12)
    def physical(x):
        return np.array([np.exp(np.log(limits[j,0])+x[j]*np.log(limits[j,1]/limits[j,0])) if j!=2 else limits[j,0]+x[j]*(limits[j,1]-limits[j,0]) for j in range(4)])
    configs=[];ids=[];unique={}
    for tr in range(8):
        row=[]
        for j in range(5):
            values=physical(points[tr,j]);config=dict(zip(names,values.tolist()))
            key=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest();unique[key]=config;row.append(key)
        ids.append(row)
    policy={'purpose':'S20 development Morris; not final GP test','names':names,'limits':limits.tolist(),'scales':['log','log','linear','log'],
        'trajectories':8,'levels':6,'delta':.6,'seed':83031,'points':points.tolist(),'orders':orders.tolist(),'query_ids':ids,
        'dx_m':args.dx,'horizon_s':args.horizon,'response':'target max event depth minus0.1m; peak depth and area also retained',
        'unknown_policy':'failure/out_of_domain is missing, never negative label; report incomplete trajectories',
        'workers':args.workers,'backend':args.backend,'scheme_order':args.order,'data_scope':'D_ONLY fixed-bed conditional propagation'}
    write_json(out/'MORRIS_POLICY.json',policy)
    physics=out/'physics';physics.mkdir();runs={}
    def query(key,config):
        folder=physics/key[:16];argv=[sys.executable,'-B',str(Path(__file__).resolve()),'pilot','--out',str(folder),
            '--inputs',str(args.inputs),'--dx',str(args.dx),'--order',str(args.order),'--backend',args.backend,
            '--horizon',str(args.horizon),'--cfl',str(args.cfl),'--dry',str(args.dry),'--maxdt',str(args.maxdt),
            '--volume',repr(config['volume_m3']),'--duration',repr(config['duration_s']),'--mu',repr(config['mu']),'--xi',repr(config['xi_m_s2']),
            '--speed',str(args.speed),'--terrain-name',args.terrain_name]
        start=time.perf_counter()
        with (physics/(key[:16]+'.stdout')).open('w') as stdout,(physics/(key[:16]+'.stderr')).open('w') as stderr:
            process=subprocess.run(argv,stdout=stdout,stderr=stderr,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        row={'query_id':key,'config':config,'argv':argv,'exit_code':process.returncode,'total_wall_s':time.perf_counter()-start,'folder':str(folder)}
        result=folder/f'pilot_{args.dx:g}.json'
        if process.returncode==0 and result.exists():
            d=json.loads(result.read_text());valid=d['whole_window_valid']
            row.update({'valid':valid,'status':d['status'],'margin_m':d['target_max_h']-.1 if valid else None,'peak_h_m':d['peak_h_anywhere_m'],'area_m2':d['inundated_area_m2'],
                'actual_volume_in_m3':d['total_inflow_m3'],'declared_volume_m3':config['volume_m3'],'solver_wall_s':d['wall_s'],'mass_error_rel':d['mass_error_rel'],'momentum_error_rel':d['momentum_error_rel'],
                'receipt':fp(folder/'RUN_RECEIPT.json')})
        else:row.update({'valid':False,'status':'numerical_failure','margin_m':None})
        return row
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(query,key,config):key for key,config in unique.items()}
        for done in as_completed(futures):
            row=done.result();runs[row['query_id']]=row
            print(json.dumps({'completed':len(runs),'unique_queries':len(unique),'query':row['query_id'][:16],'status':row['status'],'margin_m':row['margin_m'],'wall_s':row['total_wall_s']}),flush=True)
    write_json(out/'MORRIS_QUERIES.json',list(runs.values()))
    response=np.array([[runs[key]['margin_m'] if runs[key]['valid'] else np.nan for key in row] for row in ids])
    complete=np.isfinite(response).all(axis=1)
    if complete.any():
        ee,summary=morris_effects(points[complete],response[complete]);np.savez(out/'MORRIS_EFFECTS.npz',points=points,orders=orders,response=response,complete=complete,ee=ee)
        write_json(out/'MORRIS_SUMMARY.json',{'names':names,**summary,'complete_trajectories':int(complete.sum()),'total_trajectories':8,'role':'exploratory domain-dependent screening, not Sobol','independent_review':'pending'})
    else:write_json(out/'MORRIS_SUMMARY.json',{'complete_trajectories':0,'status':'pending','reason':'no complete valid trajectory; no imputation'})


class SingleGP:

    def fit(self,x,y):
        self.x=np.asarray(x);self.mean=float(np.mean(y));self.scale=max(float(np.std(y)),1e-3);self.y=(np.asarray(y)-self.mean)/self.scale
        if self.x.ndim!=2 or len(self.x)!=len(self.y) or not np.isfinite(self.x).all() or not np.isfinite(self.y).all():raise ValueError('finite matching single GP data required')
        self.noise=1e-8;dim=self.x.shape[1]
        def matrix(theta):return matern52(self.x,self.x,np.exp(theta[0]),np.exp(theta[1:]))+np.eye(len(y))*self.noise
        def loss(theta):
            try:
                c=cho_factor(matrix(theta),lower=True);return .5*self.y@cho_solve(c,self.y)+np.log(np.diag(c[0])).sum()+.5*len(y)*np.log(2*np.pi)
            except np.linalg.LinAlgError:return 1e30
        attempts=[minimize(loss,np.r_[0,np.full(dim,np.log(l))],method='L-BFGS-B',bounds=[(-3,3)]+[(-3,1)]*dim,options={'maxiter':100,'ftol':1e-8}) for l in (.15,.4,.8)]
        valid=[]
        for a in attempts:
            if not np.isfinite(a.x).all() or not np.isfinite(a.fun) or a.fun>=1e29:continue
            try:
                K=matrix(a.x)
                if np.isfinite(K).all():cho_factor(K,lower=True);valid.append(a)
            except np.linalg.LinAlgError:pass
        if not valid:raise ArithmeticError('no finite SPD single GP optimizer iterate')
        converged=[a for a in valid if a.success];best=min(converged or valid,key=lambda a:a.fun);self.theta=best.x;self.cf=cho_factor(matrix(best.x),lower=True);self.alpha=cho_solve(self.cf,self.y)
        self.fit_record={'theta':self.theta.tolist(),'mean':self.mean,'scale':self.scale,'jitter_normalized_variance':self.noise,
            'NLL':float(best.fun),'covariance_condition_number':float(np.linalg.cond(matrix(best.x))),
            'fit_status':'converged_MLE' if converged else 'finite_nonconverged_fallback',
            'fallback_policy':'prefer finite SPD converged starts; otherwise flag lowest finite SPD iterate',
            'optimizer_attempts':[{'success':bool(a.success),'fun':float(a.fun) if np.isfinite(a.fun) else None,'message':str(a.message)} for a in attempts]}
        return self
    def predict(self,x):
        k=matern52(np.asarray(x),self.x,np.exp(self.theta[0]),np.exp(self.theta[1:]))
        var=np.exp(2*self.theta[0])-np.sum(k*cho_solve(self.cf,k.T).T,axis=1)
        if var.min()<-1e-7:raise ArithmeticError('single GP negative variance')
        return self.mean+self.scale*k@self.alpha,self.scale*np.sqrt(np.maximum(0,var))


def gp_response(depth,harr=.1):
    h=np.asarray(depth,dtype=float)
    if not np.isfinite(h).all() or np.any(h<0) or not np.isfinite(harr) or harr<=0:raise ValueError('finite nonnegative depth and positive threshold required')
    return np.log1p((h-harr)/(2*harr))


def solver_identity():

    import inspect,ast,hashlib
    names=['hll_pair','minmod3','muscl_faces','exact_wall_pressure','triangle_cdf','pulse_volume','euler_step','step','energy','evolve','sample_terrain','pilot']
    functions={name:hashlib.sha256(ast.dump(ast.parse(inspect.getsource(globals()[name])),include_attributes=False).encode()).hexdigest() for name in names}
    binding={'G':G,'functions':functions,'C_source_sha256':hashlib.sha256(SOURCE.encode()).hexdigest()}
    return {'sha256':hashlib.sha256(json.dumps(binding,sort_keys=True).encode()).hexdigest(),'details':binding}


def algorithm_identity():
    import inspect,ast,hashlib
    names=['matern52','ARGP','SingleGP','gp_response','config_from_plane','execute_physics','PhysicsRegistry','save_gp_state','method_predictions','validate_execution','validate_prediction_models','run_formal_training','run_formal_test','seal_cached_seed','load_fixed_gp','cache_logical_id','FollowupRegistry','import_seed137','run_followup_validation','evaluate_followup','seal_followup_subset']
    details={name:hashlib.sha256(ast.dump(ast.parse(inspect.getsource(globals()[name])),include_attributes=False).encode()).hexdigest() for name in names}
    return {'sha256':hashlib.sha256(json.dumps(details,sort_keys=True).encode()).hexdigest(),'details':details}


def validate_execution(policy,args,policy_path):
    expected_methods={'A_single_fill','B_MF_fill','Bm_MF_count_matched','C_MF_active'}
    if set(policy['methods'])!=expected_methods or len(policy['methods'])!=4 or policy['seeds'] not in ([137],[137,271,503]):raise ValueError('exact declared method/seed set required')
    if policy['solver_identity']['sha256']!=solver_identity()['sha256'] or policy['algorithm_identity']['sha256']!=algorithm_identity()['sha256']:raise ValueError('frozen solver/algorithm changed')
    for name,expected in policy['execution'].items():
        if getattr(args,name)!=expected:raise ValueError('frozen numerical setting changed: '+name)
    for asset in policy['input_assets']:
        p=Path(args.inputs)/asset['logical_id']
        if not p.is_file() or fp(p)['sha256']!=asset['sha256']:raise ValueError('frozen input asset changed: '+asset['logical_id'])
    if policy['h_arr_m']!=.1 or policy['horizon_s']!=args.horizon:raise ValueError('arrival/window semantics changed')
    expected_geometry={'entry_x_m':350640,'entry_y_interval_m':[3135025,3135055],'target_axis':'x','target_coordinate_m':348720,'target_tangent_interval_m':[3134880,3135600]}
    if policy['geometry']!=expected_geometry or policy['low_dx_m']!=40 or policy['high_dx_m']!=30:raise ValueError('frozen geometry/fidelity semantics changed')
    for dependency in policy['dependencies']:bound_file(policy_path,dependency)
    td=json.loads(bound_file(policy_path,policy['test_design']).read_text())
    if policy['test_global_count']!=16 or policy['test_boundary_count']!=516:raise ValueError('new16 global and516 static potential points required')
    global_points=np.asarray(td['global']);critical=np.asarray(td['development_predicted_boundary'])
    expected=np.array([(z,j/128) for z in (.125,.375,.625,.875) for j in range(129)])
    if global_points.shape!=(16,2) or not np.isfinite(global_points).all() or np.any((global_points<0)|(global_points>1)):raise ValueError('invalid global coordinates')
    if not np.array_equal(critical,expected):raise ValueError('static critical tree changed')
    return True


def validate_prediction_models(methods,registry_rows,test,policy,training_root):
    expected={(int(seed),tag) for seed in policy['seeds'] for tag in policy['methods']}
    actual=[(int(m['seed']),m['method']) for m in methods]
    if len(actual)!=len(set(actual)) or set(actual)!=expected:raise ValueError('duplicate/missing/unexpected seed-method ID')
    records={r['key']:r for r in registry_rows};by_id={(int(m['seed']),m['method']):m for m in methods}
    for meta in methods:
        logical=Path(meta['folder']);folder=Path(training_root)/logical
        if logical.is_absolute() or not folder.resolve().is_relative_to(Path(training_root).resolve()):raise ValueError('method folder must be relative inside training root')
        state=json.loads((folder/'GP_STATE.json').read_text())
        if state.get('seed')!=meta['seed'] or state.get('method')!=meta['method']:raise ValueError('model seed/method mismatch')
        for kind,label in [('low','L'),('high','H')]:
            keys=meta[kind+'_keys']
            if state[kind+'_keys']!=keys or len(keys)!=len(set(keys)):raise ValueError('model training ID mismatch/duplicate')
            if any(k not in records or not records[k]['valid'] or records[k]['level']!=label for k in keys):raise ValueError('invalid or wrong-layer training ID')
            if keys:
                x=np.array([records[k]['x'] for k in keys]);y=np.array([records[k]['response_g'] for k in keys])
                if not np.array_equal(np.array(state['x'+label]),x) or not np.array_equal(np.array(state['y'+label]),y):raise ValueError('model training data not bound to physical registry')
        fit=state['fit'];theta=np.array(fit['theta']);nll=fit.get('negative_log_likelihood',fit.get('NLL'))
        if not np.isfinite(theta).all() or not np.isfinite(nll) or fit['fit_status'] not in ('converged_MLE','finite_nonconverged_fallback'):raise ValueError('invalid frozen GP fit')
        f=np.load(folder/'TEST_PREDICTIONS_FROZEN.npz')
        if f['x'].shape!=test.shape or not np.array_equal(f['x'],test):raise ValueError('prediction coordinate order mismatch')
        if f['mean_g'].shape!=(len(test),) or f['std_g'].shape!=(len(test),) or not np.isfinite(f['mean_g']).all() or not np.isfinite(f['std_g']).all() or np.any(f['std_g']<0):raise ValueError('invalid prediction shape/finite/std')
    for seed in policy['seeds']:
        A=by_id[(seed,'A_single_fill')];B=by_id[(seed,'B_MF_fill')];Bm=by_id[(seed,'Bm_MF_count_matched')];C=by_id[(seed,'C_MF_active')]
        if A['high_keys']!=B['high_keys'] or B['high_keys']!=Bm['high_keys'] or len(A['high_keys'])!=16 or len(C['high_keys'])!=16:raise ValueError('HF matching budget/design mismatch')
        if len(A['low_keys']) or len(B['low_keys'])!=24 or len(Bm['low_keys'])!=len(C['low_keys']) or not 24<=len(C['low_keys'])<=32:raise ValueError('LF count matching mismatch')
    return True


def config_from_plane(point,policy):
    config=dict(policy['fixed_parameters'])
    for x,item in zip(point,policy['axes']):
        lo,hi=item['bounds'];config[item['name']]=float(np.exp(np.log(lo)+x*np.log(hi/lo))) if item['scale']=='log' else float(lo+x*(hi-lo))
    return config


def execute_physics(folder,config,args,dx=None):

    import subprocess,sys,os
    folder=Path(folder);folder.parent.mkdir(parents=True,exist_ok=True)
    dx=float(dx if dx is not None else config.get('dx_m',args.dx))
    inputs=Path(config.get('input_dir',args.inputs))
    argv=[sys.executable,'-B',str(Path(__file__).resolve()),'pilot','--out',str(folder),'--inputs',str(inputs),
        '--dx',repr(dx),'--order',str(args.order),'--backend',args.backend,'--horizon',str(config.get('horizon_s',args.horizon)),
        '--cfl',str(config.get('cfl',args.cfl)),'--dry',str(config.get('dry_m',args.dry)),'--maxdt',str(args.maxdt),
        '--volume',repr(config['volume_m3']),'--duration',repr(config['duration_s']),'--mu',repr(config['mu']),'--xi',repr(config['xi_m_s2']),
        '--entry-width',str(config.get('entry_width_m',args.entry_width)),'--speed',str(config.get('speed_m_s',args.speed)),'--shape',config.get('shape','triangle'),'--terrain-name',config.get('terrain_name','copernicus')]
    if config.get('full',False):argv+=['--full']
    started=time.perf_counter()
    with folder.with_suffix('.stdout').open('w') as stdout,folder.with_suffix('.stderr').open('w') as stderr:
        proc=subprocess.run(argv,stdout=stdout,stderr=stderr,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    row={'config':config,'dx_m':dx,'argv':argv,'exit_code':proc.returncode,'total_wall_s':time.perf_counter()-started,'folder':str(folder)}
    result=folder/f'pilot_{dx:g}.json'
    if proc.returncode or not result.exists():
        row.update({'valid':False,'status':'numerical_failure'});return row
    d=json.loads(result.read_text());arr=d['threshold_arrivals']['0.1'];valid=d['whole_window_valid']
    row.update({'valid':valid,'status':d['status'],'target_max_h_m':d['target_max_h'],'response_g':float(gp_response(d['target_max_h'])) if valid else None,
        'raw_margin_m':d['target_max_h']-.1 if valid else None,'response_g_unit':'dimensionless_log_margin',
        'time_s':float(np.mean(arr['time_bracket_s'])) if arr['status']=='arrived' else None,'time_bracket_s':arr['time_bracket_s'],
        'solver_wall_s':d['wall_s'],'mass_error_rel':d['mass_error_rel'],'momentum_error_rel':d['momentum_error_rel'],
        'receipt':fp(folder/'RUN_RECEIPT.json'),'result':fp(result)})
    receipt=json.loads((folder/'RUN_RECEIPT.json').read_text());row['compile_and_load_s']=float((receipt.get('acceleration') or {}).get('compile_and_load_s',0.))
    return row


def run_near_boundary(out,args):
    config={'volume_m3':4e5,'duration_s':960.,'mu':.03,'xi_m_s2':500.,'speed_m_s':30.}
    history=[]
    def low(mu,index):
        c={**config,'mu':float(mu)};row=execute_physics(out/f'coarse_{index:02d}',c,args,40.)
        history.append(row);write_json(out/'COARSE_HISTORY.json',history)
        if not row['valid']:raise RuntimeError('invalid coarse bracket query; not a negative sample')
        return row
    a,b=.03,.15;left=low(a,0);right=low(b,1)
    if not (left['response_g']>=0 and right['response_g']<0):
        write_json(out/'NEAR_POINT_FROZEN.json',{'status':'no_bracket','history':history,'no_fine_queries':True});return
    for i in range(8):
        mid=.5*(a+b);row=low(mid,i+2)
        if row['response_g']>=0:a=mid;left=row
        else:b=mid;right=row
    c={**config,'mu':.5*(a+b)};frozen={'config':c,'mu_bracket':[a,b],'left':left,'right':right,'selection':'8 fixed40m bisections, no fine label accessed','fine_status':'not_run'}
    write_json(out/'NEAR_POINT_FROZEN.json',frozen);binding=fp(out/'NEAR_POINT_FROZEN.json');print(json.dumps({'frozen':binding,'config':c}),flush=True)
    fine30=execute_physics(out/'fine30',c,args,30.);fine15=execute_physics(out/'fine15',c,args,15.)
    write_json(out/'NEAR_COMPARISON.json',{'frozen_point':binding,'coarse_bracket':[left,right],'fine30':fine30,'fine15':fine15,
        'scope':'one development near-LF-boundary; numerical ambiguity retained','independent_review':'pending'})


def run_perturbations(out,args):
    from concurrent.futures import ThreadPoolExecutor,as_completed
    base={'volume_m3':args.volume,'duration_s':args.duration,'mu':args.mu,'xi_m_s2':args.xi,'speed_m_s':args.speed}
    cases={'base':base,'T720':{**base,'duration_s':720.},'T1200':{**base,'duration_s':1200.},'xi250':{**base,'xi_m_s2':250.},'xi1000':{**base,'xi_m_s2':1000.}}
    write_json(out/'PERTURBATION_DESIGN.json',{'cases':cases,'dx_m':30.,'role':'development-only confirmation of fixed parameters','test_truth_accessed':False})
    result={}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(execute_physics,out/name,c,args,30.):name for name,c in cases.items()}
        for f in as_completed(futures):
            name=futures[f];result[name]=f.result();print(json.dumps({'case':name,'status':result[name]['status'],'g':result[name].get('response_g')}),flush=True)
    write_json(out/'PERTURBATION_RESULTS.json',result)


def run_global_cases(out,args):

    from concurrent.futures import ThreadPoolExecutor,as_completed
    common={'volume_m3':8e5,'duration_s':960.,'mu':.05,'xi_m_s2':500.,'speed_m_s':30.,'full':True,'horizon_s':3600.,'input_dir':str(args.full_inputs)}
    cases={}
    for terrain_name in ('copernicus','aster'):
        for shape,duration in [('triangle',960.),('wide',1200.),('double',960.)]:
            cases[terrain_name+'_'+shape+'_40']={**common,'terrain_name':terrain_name,'shape':'triangle' if shape=='wide' else shape,'duration_s':duration,'dx_m':40.}
    cases['copernicus_triangle_30']={**common,'terrain_name':'copernicus','shape':'triangle','dx_m':30.}
    cases['copernicus_high_mobility_40']={**common,'volume_m3':8e5,'mu':.03,'xi_m_s2':1000.,'terrain_name':'copernicus','shape':'triangle','dx_m':40.}
    cases['copernicus_low_mobility_40']={**common,'volume_m3':2e5,'mu':.15,'xi_m_s2':250.,'terrain_name':'copernicus','shape':'triangle','dx_m':40.}
    write_json(out/'GLOBAL_CASE_POLICY.json',{'cases':cases,'same_volume_m3':8e5,'global_horizon_s':3600,
        'two_pulses':'two adjacent V/2,T/2 triangles inside [0,T], same total support/peak as single triangle',
        'wide_pulse':'sameV but longerT; pressure/energy input not assumed identical','scope':'full ~14.188km conditional reach, not local GP accuracy or event reconstruction',
        'unknown_observed_entry_time':None,'independent_review':'pending'})
    result={}
    with ThreadPoolExecutor(max_workers=min(args.workers,2)) as pool:
        futures={pool.submit(execute_physics,out/name,c,args,c['dx_m']):name for name,c in cases.items()}
        for done in as_completed(futures):
            name=futures[done];result[name]=done.result();print(json.dumps({'full_case':name,'status':result[name]['status'],'target_max_h_m':result[name].get('target_max_h_m'),'wall_s':result[name]['total_wall_s']}),flush=True)
            write_json(out/'GLOBAL_RESULTS.json',result)
    table=[]
    for name,row in result.items():
        detail=json.loads(Path(row['result']['path']).read_text()) if row['valid'] else {}
        table.append({'case_id':name,**{k:row['config'][k] for k in ['volume_m3','duration_s','mu','xi_m_s2','terrain_name','shape']},
            'H_s':3600,'dx_m':row['dx_m'],'status':row['status'],'travel_time_s':row.get('time_s'),'target_max_h_m':row.get('target_max_h_m'),
            'domain_space_time_peak_h_m':detail.get('peak_h_anywhere_m'),'ever_h_ge0p1_coverage_area_m2':detail.get('ever_threshold_coverage_area_m2'),
            'actual_inflow_m3':detail.get('total_inflow_m3'),'mass_error_rel':row.get('mass_error_rel'),'query_wall_s':row['total_wall_s']})
    import csv
    with (out/'GLOBAL_TABLE.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)


def along_axis_baseline(inputs,full,mu,xi,h_vertical=2.,u0=30.,horizon=3600.):

    from scipy.integrate import quad
    from scipy.optimize import brentq
    from scipy.ndimage import map_coordinates
    coordinates=np.loadtxt(Path(inputs)/'mapped_path_utm.csv',delimiter=',',skiprows=1)
    terrain=np.load(Path(inputs)/'copernicus_support.npz');bounds=terrain['bounds'];res=float(terrain['resolution']);bed=terrain['z']
    distance=np.r_[0,np.cumsum(np.linalg.norm(np.diff(coordinates,axis=0),axis=1))]
    def crossing(axis,value,start=0):
        for j in range(start,len(coordinates)-1):
            a,b=coordinates[j,axis],coordinates[j+1,axis]
            if a!=b and (a-value)*(b-value)<=0:
                f=(value-a)/(b-a);return j,float(distance[j]+f*(distance[j+1]-distance[j]))
        raise ValueError('mapped axis does not intersect the declared boundary/target')
    i0,s0=crossing(0,350640);it,st=crossing(1,3129103.736725446,i0) if full else crossing(0,348720,i0)
    s=np.r_[s0,distance[(distance>s0)&(distance<st)],st]-s0
    x=np.interp(s+s0,distance,coordinates[:,0]);y=np.interp(s+s0,distance,coordinates[:,1])
    z=map_coordinates(bed,[(bounds[3]-y)/res-.5,(x-bounds[0])/res-.5],order=1,mode='nearest')
    energy=.5*u0*u0;t=0.;position=0.;stopped=False;trajectory=[[0.,0.,u0]]
    for j,ds in enumerate(np.diff(s)):
        slope=(z[j]-z[j+1])/ds;a=G*(slope-mu);b=2*G/(xi*h_vertical);eq=a/b
        if energy<=0 and a<=0:stopped=True;break
        e0=energy
        def e(l):return eq+(e0-eq)*np.exp(-b*l)
        reach=float(ds)
        if e(reach)<=0 and eq<0:reach=float(-np.log(-eq/(e0-eq))/b);stopped=True
        def travel(l):return quad(lambda u:1/np.sqrt(max(2*e(u),1e-20)),0,l,epsabs=1e-7,epsrel=1e-7,limit=100)[0]
        elapsed=travel(reach)
        if t+elapsed>horizon:
            reach=brentq(lambda l:travel(l)-(horizon-t),0,reach);position+=reach;t=horizon;energy=max(0,float(e(reach)));trajectory.append([t,position,math.sqrt(2*energy)]);break
        position+=reach;t+=elapsed;energy=max(0,float(e(reach)));trajectory.append([t,position,math.sqrt(2*energy)])
        if stopped:break
    arrived=position>=s[-1]-1e-6
    return {'status':'arrived' if arrived else 'valid_no_arrival_by_horizon','travel_time_s':t if arrived else None,'axis_position_at_cutoff_or_stop_m':position,
        'target_axis_distance_m':float(s[-1]),'baseline_stopped':stopped,'trajectory':trajectory,'h_vertical_closure_m':h_vertical,
        'scope':'horizontal along-axis energy ODE, constant assumed thickness, no mass hydrograph/lateral inundation; not low fidelity or event validation'}


def run_synthetic_inverse(out,args):

    from concurrent.futures import ThreadPoolExecutor,as_completed
    truth={'volume_m3':7e5,'duration_s':960.,'mu':.04,'xi_m_s2':500.,'speed_m_s':30.,'cfl':.175}
    vv=np.linspace(6e5,8e5,4);mm=np.linspace(.03,.05,4)
    configurations=[{'volume_m3':float(v),'duration_s':960.,'mu':float(mu),'xi_m_s2':500.,'speed_m_s':30.} for v in vv for mu in mm]
    policy={'truth':truth,'observation_generator_dx_m':15.,'candidate_model_dx_m':30.,'candidate_library':configurations,
        'arrival_error_scale_s':60.,'target_depth_error_scale_m':.75,'candidate_parameter_grid_per_axis':41,
        'calibration_outputs':['arrival_time','target_peak_depth'],'held_out_output':'target_peak_signed_discharge',
        'scope':'synthetic observation with grid/time mismatch, not local calibration or natural error distribution',
        'representatives':'min arrival error, min depth error, min max(normalized errors), unique points returned to PDE'}
    write_json(out/'INVERSE_POLICY.json',policy)

    observed=execute_physics(out/'synthetic_truth_15',truth,args,15.)
    if not observed['valid'] or observed['time_s'] is None:raise RuntimeError('registered synthetic truth did not arrive; keep failure, do not choose another truth')
    rows=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(execute_physics,out/f'library_{j:02d}',c,args,30.):j for j,c in enumerate(configurations)}
        by_index={}
        for done in as_completed(futures):
            j=futures[done];by_index[j]=done.result();print(json.dumps({'inverse_training':j,'status':by_index[j]['status']}),flush=True)
    rows=[by_index[j] for j in range(len(configurations))]
    write_json(out/'INVERSE_PHYSICAL_LIBRARY.json',{'observed':observed,'training':rows})
    if any(not r['valid'] or r['time_s'] is None for r in rows):raise RuntimeError('inverse predeclared library contains nonarrival/unknown; do not impute finite time')
    x=np.array([[(r['config']['volume_m3']-6e5)/2e5,(r['config']['mu']-.03)/.02] for r in rows])
    times=np.array([r['time_s'] for r in rows]);depths=np.array([r['target_max_h_m'] for r in rows])
    gp_t=SingleGP().fit(x,times);gp_h=SingleGP().fit(x,depths)
    grid=np.array([(a,b) for a in np.linspace(0,1,41) for b in np.linspace(0,1,41)])
    pred_t,_=gp_t.predict(grid);pred_h,_=gp_h.predict(grid)
    e1=abs(pred_t-observed['time_s'])/60.;e2=abs(pred_h-observed['target_max_h_m'])/.75
    order=np.argsort(e1,kind='stable');front=[];best=np.inf
    for j in order:
        if e2[j]<best:front.append(int(j));best=e2[j]
    representatives=list(dict.fromkeys([int(np.argmin(e1)),int(np.argmin(e2)),int(np.argmin(np.maximum(e1,e2)))]))
    checked=[];observed_detail=json.loads(Path(observed['result']['path']).read_text())
    for i in representatives:
        c={'volume_m3':float(6e5+2e5*grid[i,0]),'duration_s':960.,'mu':float(.03+.02*grid[i,1]),'xi_m_s2':500.,'speed_m_s':30.}
        actual=execute_physics(out/f'pareto_recheck_{i:04d}',c,args,30.)
        detail=json.loads(Path(actual['result']['path']).read_text()) if actual['valid'] else {}
        checked.append({'grid_index':i,'normalized_parameter':grid[i].tolist(),'surrogate_errors':[float(e1[i]),float(e2[i])],'physical':actual,
            'held_out_Q_error_m3_s':detail.get('target_Q_downstream_max_m3_s',np.nan)-observed_detail['target_Q_downstream_max_m3_s'] if detail else None,
            'truth_parameter_relative_errors':{'volume':(c['volume_m3']/truth['volume_m3']-1),'mu':(c['mu']/truth['mu']-1)}})
    np.savez_compressed(out/'INVERSE_SEARCH.npz',training_x=x,training_time_s=times,training_depth_m=depths,grid=grid,predicted_time_s=pred_t,predicted_depth_m=pred_h,
        normalized_arrival_error=e1,normalized_depth_error=e2,pareto_indices=front,single_compatible=e1<=1,multi_compatible=(e1<=1)&(e2<=1))

    slope=.10;h=2.;mu0=.05;xi0=500.;terminal=math.sqrt(xi0*h*(slope-mu0));tobs=10.
    transient=terminal*math.tanh(G*(slope-mu0)*tobs/terminal)
    family=[{'mu':mu,'xi':terminal**2/(h*(slope-mu)),'terminal_m_s':terminal,
        'u_at10s_m_s':terminal*math.tanh(G*(slope-mu)*tobs/terminal)} for mu in (.02,.035,.05,.065,.075)]
    recovered_mu=slope-terminal/(G*tobs)*math.atanh(transient/terminal);recovered_xi=terminal**2/(h*(slope-recovered_mu))
    write_json(out/'INVERSE_RESULTS.json',{'policy':fp(out/'INVERSE_POLICY.json'),'observed':observed,'arrival_GP':gp_t.fit_record,'depth_GP':gp_h.fit_record,
        'single_constraint_candidates':int((e1<=1).sum()),'two_constraint_candidates':int(((e1<=1)&(e2<=1)).sum()),
        'candidate_grid_count':len(grid),'Pareto_indices':front,'representatives_recomputed':checked,
        'nonidentifiable_terminal_velocity_family':family,'two_observable_exact_recovery':{'true_mu':mu0,'recovered_mu':recovered_mu,'true_xi':xi0,'recovered_xi':recovered_xi},
        'not_a_natural_calibration':True,'independent_review':'pending'})


def strict_q3_q4_join(q3,q4,entity_root,allow_synthetic=False):


    entity_root=Path(entity_root).resolve()
    if q3.get('expected_targets')!=q4.get('expected_targets'):raise ValueError('expected target geometries differ')
    expected=q4.get('expected_targets',{})
    if not expected:raise ValueError('expected targets required')
    def verify(record):
        p=(entity_root/record['path']).resolve()
        if not p.is_relative_to(entity_root) or not p.is_file() or fp(p)['sha256']!=record['sha256']:raise ValueError('entity missing, outside root or fingerprint mismatch')
        return p
    entries={}
    for entry in q3.get('entries',[]):
        sid=entry['scenario_id']
        if sid in entries:raise ValueError('duplicate Q3 scenario')
        if entry.get('entry_stage')!='propagation_entry':raise ValueError('local first hit is not automatically propagation entry')
        bridge=entry.get('bridge_evidence')
        if not bridge or not entry.get('bridge_evidence_id'):raise ValueError('bound stage-bridge evidence required, not only a label')
        b=json.loads(verify(bridge).read_text())
        if b.get('source_kind')=='synthetic_assumption' and not allow_synthetic:raise ValueError('synthetic stage assumption cannot enter a natural join')
        if b.get('source_kind') not in ('synthetic_assumption','independently_reviewed_stage_bridge'):raise ValueError('unreviewed stage bridge')
        if entry.get('status') not in ('entry_reached','not_entered','unknown'):raise ValueError('invalid Q3 entry status')
        entries[sid]=entry
    responses={};identities=['scenario_hash','joint_input_hash','entry_id','entry_geometry_hash','entry_state_hash','time_origin_id']
    entity_fields={'parameters':'parameter_hash','terrain':'terrain_hash','code':'code_hash','entry_geometry':'entry_geometry_hash','entry_state':'entry_state_hash',
        'target_geometry':'target_geometry_hash','criterion':'arrival_criterion_hash','distance':'distance_definition_hash'}
    for response in q4.get('arrivals',[]):
        key=(response['scenario_id'],response['target_id'])
        if key in responses or key[0] not in entries or key[1] not in expected:raise ValueError('duplicate or unexpected Q4 scenario-target')
        entry=entries[key[0]]
        if any(entry.get(k) is None or entry.get(k)!=response.get(k) for k in identities):raise ValueError('joint member, entry/state or clock mismatch')
        if response['target_geometry_hash']!=expected[key[1]] or response.get('time_unit')!='s' or response.get('distance_unit')!='m':raise ValueError('target geometry or SI-unit mismatch')
        entities=response.get('entities',{})
        for kind,field in entity_fields.items():
            if kind not in entities:raise ValueError('missing bound entity: '+kind)
            verify(entities[kind])
            if response.get(field)!=entities[kind]['sha256']:raise ValueError('declared entity identity differs: '+kind)
        criterion=json.loads(verify(entities['criterion']).read_text());distance=json.loads(verify(entities['distance']).read_text())
        if criterion.get('definition')!=response.get('arrival_criterion') or distance.get('definition')!=response.get('distance_definition'):raise ValueError('distance/arrival semantics mismatch')
        status=response.get('target_arrival_status');travel=response.get('travel_time_s');valid_until=response.get('valid_until_s')
        horizon=response.get('propagation_horizon_s')
        if response.get('horizon_s',horizon)!=horizon:raise ValueError('conflicting horizon aliases')
        if status=='arrived':
            if entry['status']=='not_entered':raise ValueError('not-entered scenario cannot have a propagated arrival')
            if travel is None or not math.isfinite(travel) or travel<0 or response.get('prefix_review_status')!='pass' or valid_until is None or not math.isfinite(valid_until) or valid_until<travel:raise ValueError('arrival lacks finite valid independently reviewed prefix')
        elif status=='valid_no_arrival_by_horizon':
            if travel is not None or horizon is None or not math.isfinite(horizon) or horizon<=0 or response.get('whole_window_valid') is not True or valid_until is None or not math.isfinite(valid_until) or valid_until<horizon:raise ValueError('nonarrival requires complete valid window')
        elif status not in ('unknown','not_run') or travel is not None:raise ValueError('invalid unknown/nonarrival time')
        responses[key]=response
    joined=[]
    for sid,entry in entries.items():
        for target in expected:
            response=responses.get((sid,target));row={'scenario_id':sid,'target_id':target,'absolute_arrival_s':None,'effective_lead_s':None,'travel_time_s':None,'target_arrival_status':'unknown'}
            if response is None:row['reason']='missing expected scenario-target response';joined.append(row);continue
            row.update(travel_time_s=response['travel_time_s'],target_arrival_status=response['target_arrival_status'],run_termination_status=response.get('run_termination_status'),
                whole_window_valid=response.get('whole_window_valid'),valid_until_s=response.get('valid_until_s'),reason='relative propagation response only')
            if response['target_arrival_status']=='arrived' and entry['status']=='entry_reached' and entry.get('T_in_s') is not None:
                if not math.isfinite(entry['T_in_s']):raise ValueError('nonfinite entry time')
                row['absolute_arrival_s']=entry['T_in_s']+response['travel_time_s']
                pub,action=entry.get('t_pub_s'),entry.get('d_act_s')
                if pub is not None and action is not None:
                    if not math.isfinite(pub) or not math.isfinite(action) or action<0:raise ValueError('invalid publication/action delay')
                    row['effective_lead_s']=row['absolute_arrival_s']-pub-action
            joined.append(row)
    return {'joined':joined,'scientific_bridge_truth':'must be independently established; entity consistency alone does not establish it'}


def freeze_formal_policy(out,args):

    import shutil
    from scipy.stats import qmc
    prior=Path(args.previous_design);policy=json.loads((prior/'FORMAL_POLICY.json').read_text());designs=json.loads((prior/'TRAINING_DESIGNS.json').read_text())
    if set(designs)!={'137','271','503'}:raise ValueError('original three-seed designs required')
    (out/'TRAINING_DESIGNS.json').write_bytes((prior/'TRAINING_DESIGNS.json').read_bytes())
    old_path=prior/'TEST_DESIGN_ORIGINAL96.json';original=json.loads(old_path.read_text());forbidden=[*original['global'],*original['development_predicted_boundary']]
    for d in designs.values():forbidden.extend(d['low']);forbidden.extend(d['active_candidates']);forbidden.extend(d['matched_passive_low_pool'])
    pool=qmc.Sobol(2,scramble=True,seed=93071).random_base2(6);f=np.asarray(forbidden);keep=[i for i,x in enumerate(pool) if np.min(np.linalg.norm(x-f,axis=1))>1e-12][:16]
    if len(keep)!=16:raise ValueError('not enough independent global coordinates')
    critical=np.array([(z,j/128) for z in (.125,.375,.625,.875) for j in range(129)])
    if np.min(np.linalg.norm(critical[:,None,:]-f[None,:,:],axis=-1))<=1e-12:raise ValueError('critical tree coincides with prior training/test coordinate')
    write_json(out/'TEST_DESIGN.json',{'global':pool[keep].tolist(),'development_predicted_boundary':critical.tolist(),'global_seed':93071,'kept_pool_indices':keep,'critical_label':'potential static tree, only selected HF queried; not development-predicted band','historical_test_coordinates':fp(old_path),'labels':'not_evaluated'})
    policy.update({'artifact_id':'Q4-FOLLOWUP-v1','seeds':[137,271,503],'test_global_count':16,'test_boundary_count':516,'test_global_seed':93071,
        'test_design':fp(out/'TEST_DESIGN.json'),'training_design':fp(out/'TRAINING_DESIGNS.json'),'solver_identity':solver_identity(),'algorithm_identity':algorithm_identity(),
        'design_generator':fp(__file__),'parent_policy':fp(prior/'FORMAL_POLICY.json'),'test_access':'12 models and all532 static predictions sealed before any new validation truth',
        'boundary_definition':'4 fixed zV*129 mu nodes; actual slices zV3/8,7/8; initial j0,64,128; leftmost adjacent opposite classes; at most4 integer bisections',
        'critical_max_H':14,'critical_max_L':4,'old24_role':'historical only, never imported as new test','independent_review':'pending'})
    for dep in policy['dependencies']:
        source=Path(dep['path']);shutil.copy2(source,out/source.name)
    write_json(out/'FORMAL_POLICY.json',policy)


def bound_file(policy_path,record):
    candidates=[Path(policy_path).parent/Path(record['path']).name,Path(record['path'])]
    for p in candidates:
        if p.is_file() and fp(p)['sha256']==record['sha256']:return p
    raise ValueError('bound policy dependency missing or changed: '+record['path'])


class PhysicsRegistry:
    def __init__(self,out,args,policy):
        self.out=Path(out);self.args=args;self.policy=policy;self.records={};self.requests=[]
        self.out.mkdir();self.source_hash=fp(__file__)['sha256'];self.terrain=fp(Path(args.inputs)/'copernicus_support.npz')['sha256']
    def key(self,x,level):
        import hashlib
        config=config_from_plane(x,self.policy);dx=self.policy['low_dx_m'] if level=='L' else self.policy['high_dx_m']
        body={'config':config,'dx':dx,'horizon':self.args.horizon,'dry':self.args.dry,'cfl':self.args.cfl,'backend':self.args.backend,'source':self.source_hash,'terrain':self.terrain}
        return hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest(),config,dx
    def request(self,x,level,role):
        key,config,dx=self.key(x,level);cached=key in self.records
        if not cached:
            row=execute_physics(self.out/key[:20],config,self.args,dx);row.update({'key':key,'x':np.asarray(x).tolist(),'level':level});self.records[key]=row
            write_json(self.out/'REGISTRY.json',list(self.records.values()))
        self.requests.append({'key':key,'role':role,'cache_hit':cached});return key
    def batch(self,x,level,role):
        from concurrent.futures import ThreadPoolExecutor,as_completed
        planned={};keys=[]
        for point in x:
            key,config,dx=self.key(point,level);keys.append(key)
            if key in self.records:self.requests.append({'key':key,'role':role,'cache_hit':True})
            else:planned[key]=(point,config,dx)
        with ThreadPoolExecutor(max_workers=self.args.workers) as pool:
            futures={pool.submit(execute_physics,self.out/key[:20],c,self.args,dx):key for key,(point,c,dx) in planned.items()}
            for future in as_completed(futures):
                key=futures[future];point,_,_=planned[key];row=future.result();row.update({'key':key,'x':np.asarray(point).tolist(),'level':level});self.records[key]=row
                self.requests.append({'key':key,'role':role,'cache_hit':False});write_json(self.out/'REGISTRY.json',list(self.records.values()))
                print(json.dumps({'role':role,'level':level,'new_completed':len(self.records),'status':row['status'],'g':row.get('response_g'),'wall_s':row['total_wall_s']}),flush=True)
        return keys
    def data(self,keys):
        rows=[self.records[k] for k in keys]
        if any(not r['valid'] for r in rows):raise RuntimeError('invalid physical training query; unknown is not a negative label')
        return np.array([r['x'] for r in rows]),np.array([r['response_g'] for r in rows])


def save_gp_state(path,model,keys_low,keys_high,registry,seed,method):
    state={'seed':seed,'method':method,'fit':model.fit_record,'low_keys':keys_low,'high_keys':keys_high,'scope':'conditional latent GP; no natural probability/95% claim'}
    if keys_low:
        x,y=registry.data(keys_low);state.update({'xL':x.tolist(),'yL':y.tolist()})
    x,y=registry.data(keys_high);state.update({'xH':x.tolist(),'yH':y.tolist()});write_json(path,state)


def method_predictions(folder,model,test,high_keys,registry):
    started=time.perf_counter()
    mean,std=model.predict(test)
    np.savez(folder/'TEST_PREDICTIONS_FROZEN.npz',x=test,mean_g=mean,std_g=std)

    axis=np.linspace(0,1,61);surface=np.array([(a,b) for b in axis for a in axis]);ms=[];ss=[]
    for start in range(0,len(surface),512):
        a,b=model.predict(surface[start:start+512]);ms.extend(a);ss.extend(b)
    np.savez_compressed(folder/'SURFACE_PREDICTIONS.npz',axis=axis,mean_g=np.array(ms).reshape(61,61),std_g=np.array(ss).reshape(61,61))
    write_json(folder/'PREDICTION_SCOPE.json',{'model':'arrival margin g only','time_proxy':'not_modelled; physical first-arrival brackets retained in query results','not_a_natural_probability':True})
    return time.perf_counter()-started


def run_formal_training(out,args):

    import csv
    policy_path=Path(args.policy);policy=json.loads(policy_path.read_text())
    validate_execution(policy,args,policy_path)
    design_path=bound_file(policy_path,policy['training_design']);test_path=bound_file(policy_path,policy['test_design'])
    designs=json.loads(design_path.read_text());test_design=json.loads(test_path.read_text())
    test=np.r_[np.array(test_design['global']),np.array(test_design['development_predicted_boundary'])]
    args.horizon=policy['horizon_s'];registry=FollowupRegistry(out/'physics',args,policy);registry.import_cache(args.cache_registry);methods=[]
    if args.training_dir is not None:methods=list(import_seed137(out,args,registry,test))
    write_json(out/'METHODS.json',methods)
    write_json(out/'EXECUTED_POLICY_BINDING.json',{'policy':fp(policy_path),'training_design':fp(design_path),'test_design':fp(test_path),'source':fp(__file__),'solver':solver_identity(),'test_labels':'not_accessed'})
    for seed in policy['seeds']:
        if seed==137 and args.training_dir is not None:continue
        d=designs[str(seed)];xl=np.array(d['low']);hi0=np.array(d['high_initial_indices']);hip=np.array(d['high_passive_added_indices'])

        pool=np.array(d['active_candidates']);integration=np.array(d['integration_grid'])
        matched_pool=np.array(d['matched_passive_low_pool'])
        if min(np.linalg.norm(np.r_[xl,pool,matched_pool][:,None,:]-test[None,:,:],axis=-1).ravel())<1e-12:raise ValueError('test/train or acquisition coordinate collision')
        low=registry.batch(xl,'L',f'seed{seed}_shared_L24');initial=registry.batch(xl[hi0],'H',f'seed{seed}_shared_H8')
        passive_added=registry.batch(xl[hip],'H',f'seed{seed}_passive_H8');passive=initial+passive_added
        xh,yh=registry.data(passive);xlow,ylow=registry.data(low)
        for tag in ('A_single_fill','B_MF_fill'):
            folder=out/f'seed{seed}'/tag;folder.mkdir(parents=True);started=time.perf_counter()
            model=SingleGP().fit(xh,yh) if tag.startswith('A') else ARGP().fit(xlow,ylow,xh,yh)
            fit_s=time.perf_counter()-started;save_gp_state(folder/'GP_STATE.json',model,[] if tag.startswith('A') else low,passive,registry,seed,tag)
            aux_s=method_predictions(folder,model,test,passive,registry)
            methods.append({'seed':seed,'method':tag,'folder':str(folder.relative_to(out)),'low_keys':[] if tag.startswith('A') else low,'high_keys':passive,'fit_wall_s':fit_s,'acquisition_wall_s':0.,'auxiliary_prediction_wall_s':aux_s})
            write_json(out/'METHODS.json',methods)
        folder=out/f'seed{seed}'/'C_MF_active';folder.mkdir(parents=True);lkeys=list(low);hkeys=list(initial);history=[];nlow=0;nhigh=0;iteration=0;fit_s=0.;acq_s=0.
        while nhigh<policy['active_high_added']:
            iteration+=1;xlow,ylow=registry.data(lkeys);xh,yh=registry.data(hkeys)
            started=time.perf_counter();model=ARGP().fit(xlow,ylow,xh,yh);fit_s+=time.perf_counter()-started
            save_gp_state(folder/f'action{iteration:02d}_MODEL.json',model,lkeys,hkeys,registry,seed,'C_MF_active')
            started=time.perf_counter();explore=iteration%4==0
            scores=model.acquisition(pool,integration,policy['estimated_cost_L_s'],policy['estimated_cost_H_s'],explore=explore)

            fallback=not any(np.any(v>0) for v in scores.values())
            if fallback:scores=model.acquisition(pool,integration,policy['estimated_cost_L_s'],policy['estimated_cost_H_s'],explore=True)
            eligibleL=np.min(np.linalg.norm(pool[:,None,:]-xlow[None,:,:],axis=-1),axis=1)>1e-12
            eligibleH=np.min(np.linalg.norm(pool[:,None,:]-xh[None,:,:],axis=-1),axis=1)>1e-12
            if nlow>=policy['active_low_cap']:eligibleL[:]=False
            scoreL=np.where(eligibleL,scores['L'],-np.inf);scoreH=np.where(eligibleH,scores['H'],-np.inf)
            chosen=int(np.argmax(np.r_[scoreL,scoreH]));level='L' if chosen<len(pool) else 'H';index=chosen%len(pool);point=pool[index]
            np.savez(folder/f'action{iteration:02d}_ACQUISITION.npz',candidates=pool,integration=integration,scoreL=scores['L'],scoreH=scores['H'],eligibleL=eligibleL,eligibleH=eligibleH)
            predicted_mean,predicted_std=model.predict(point[None,:]);acq_s+=time.perf_counter()-started
            key=registry.request(point,level,f'seed{seed}_active_action{iteration}');record=registry.records[key]
            if not record['valid']:raise RuntimeError('active physical query invalid; method aborted, not relabeled')
            if level=='L':lkeys.append(key);nlow+=1
            else:hkeys.append(key);nhigh+=1
            history.append({'iteration':iteration,'x':point.tolist(),'level':level,'key':key,'exploration':explore,'underflow_global_fallback':fallback,
                'predicted_H_mean':float(predicted_mean[0]),'predicted_H_std':float(predicted_std[0]),'selected_score':float(scores[level][index]),
                'cost_estimate_s':policy['estimated_cost_'+level+'_s'],'actual_total_wall_s':record['total_wall_s'],'actual_g':record['response_g']})
            write_json(folder/'ACTIVE_HISTORY.json',history);write_json(folder/'ACTIVE_TIMING.json',{'fit_wall_s':fit_s,'acquisition_wall_s':acq_s,'high_added':nhigh,'low_added':nlow});print(json.dumps({'seed':seed,'active_action':iteration,'level':level,'H_added':nhigh,'L_added':nlow}),flush=True)
        xlow,ylow=registry.data(lkeys);xh,yh=registry.data(hkeys);started=time.perf_counter();model=ARGP().fit(xlow,ylow,xh,yh);fit_s+=time.perf_counter()-started
        save_gp_state(folder/'GP_STATE.json',model,lkeys,hkeys,registry,seed,'C_MF_active');aux_s=method_predictions(folder,model,test,hkeys,registry)
        methods.append({'seed':seed,'method':'C_MF_active','folder':str(folder.relative_to(out)),'low_keys':lkeys,'high_keys':hkeys,'fit_wall_s':fit_s,'acquisition_wall_s':acq_s,'auxiliary_prediction_wall_s':aux_s})
        write_json(out/'METHODS.json',methods)


        extra=registry.batch(matched_pool[:nlow],'L',f'seed{seed}_Bm_matched_passive_L{nlow}') if nlow else []
        bm_low=low+extra;xlow,ylow=registry.data(bm_low);xh,yh=registry.data(passive)
        bm_folder=out/f'seed{seed}'/'Bm_MF_count_matched';bm_folder.mkdir();started=time.perf_counter();bm=ARGP().fit(xlow,ylow,xh,yh);bm_fit=time.perf_counter()-started
        save_gp_state(bm_folder/'GP_STATE.json',bm,bm_low,passive,registry,seed,'Bm_MF_count_matched');bm_aux=method_predictions(bm_folder,bm,test,passive,registry)
        methods.append({'seed':seed,'method':'Bm_MF_count_matched','folder':str(bm_folder.relative_to(out)),'low_keys':bm_low,'high_keys':passive,'fit_wall_s':bm_fit,'acquisition_wall_s':0.,'auxiliary_prediction_wall_s':bm_aux})
        write_json(out/'METHODS.json',methods)
    def local_fp(path):
        p=Path(path);record=fp(p);record['path']=str(p.relative_to(out));return record
    prediction_manifest=[local_fp(out/m['folder']/'TEST_PREDICTIONS_FROZEN.npz') for m in methods]
    expected=len(policy['seeds'])*len(policy['methods'])
    if len(prediction_manifest)!=expected:raise RuntimeError('all method predictions must be frozen before test truth')
    validate_prediction_models(methods,list(registry.records.values()),test,policy,out)
    if any(np.min(np.linalg.norm(np.array(r['x'])[None,:]-test,axis=1))<1e-12 for r in registry.records.values()):raise ValueError('test coordinate already present in training physical library')
    write_json(out/'METHODS.json',methods);write_json(out/'QUERY_REQUESTS.json',registry.requests)
    seal={'phase':'training_complete_test_not_run','predictions':prediction_manifest,'model_states':[local_fp(out/m['folder']/'GP_STATE.json') for m in methods],
        'policy':fp(policy_path),'test_design':fp(test_path),'source':fp(__file__),'solver':solver_identity(),'methods':local_fp(out/'METHODS.json'),
        'registry':local_fp(out/'physics/REGISTRY.json'),'requests':local_fp(out/'QUERY_REQUESTS.json'),'terrain':fp(Path(args.inputs)/'copernicus_support.npz'),
        'numerics':{'order':args.order,'backend':args.backend,'dry':args.dry,'cfl':args.cfl,'maxdt':args.maxdt},'no_test_label_coordinates_in_registry':True}
    write_json(out/'PHASE1_SEAL.json',seal);print(json.dumps({'phase1_seal':fp(out/'PHASE1_SEAL.json'),'test_queries':0}),flush=True)


def seal_cached_seed(out,args):


    import shutil,importlib.util,sys
    parent=Path(args.training_dir);policy_path=Path(args.policy);policy=json.loads(policy_path.read_text())
    validate_execution(policy,args,policy_path)
    binding=json.loads((parent/'EXECUTED_POLICY_BINDING.json').read_text())
    snapshot=parent/'SOURCE_SNAPSHOT.py'
    if fp(snapshot)['sha256']!=binding['source']['sha256']:raise ValueError('parent source snapshot changed')
    spec=importlib.util.spec_from_file_location('q4_prior_cache_source',snapshot);prior=importlib.util.module_from_spec(spec);sys.modules[spec.name]=prior;spec.loader.exec_module(prior)
    if prior.solver_identity()['sha256']!=solver_identity()['sha256'] or binding['solver']['sha256']!=solver_identity()['sha256']:raise ValueError('physical solver differs from cached labels')
    if (parent/'TEST_ACCESS_STARTED.json').exists():raise ValueError('parent test already accessed')
    td=json.loads(bound_file(policy_path,policy['test_design']).read_text());test=np.r_[td['global'],td['development_predicted_boundary']]
    all_rows=json.loads((parent/'physics/REGISTRY.json').read_text());records={r['key']:r for r in all_rows};methods=[];used=set();lineage=[]
    for tag in policy['methods']:
        relative=Path('seed137')/tag;old=parent/relative;folder=out/relative
        shutil.copytree(old,folder)
        state=json.loads((folder/'GP_STATE.json').read_text())
        oldpred=folder/'TEST_PREDICTIONS_FROZEN.npz';saved=np.load(oldpred);ox=saved['x'].copy();om=saved['mean_g'].copy();os=saved['std_g'].copy();saved.close()
        indices=[]
        for point in test:
            found=np.flatnonzero(np.all(ox==point,axis=1))
            if len(found)!=1:raise ValueError('new coordinate not uniquely in original frozen predictions')
            indices.append(int(found[0]))
        oldpred.rename(folder/'TEST_PREDICTIONS_PARENT48.npz')
        np.savez(oldpred,x=ox[indices],mean_g=om[indices],std_g=os[indices])
        keys=state['low_keys']+state['high_keys'];used.update(keys)
        methods.append({'seed':137,'method':tag,'folder':str(relative),'low_keys':state['low_keys'],'high_keys':state['high_keys'],
            'fit_wall_s':None,'acquisition_wall_s':None,'auxiliary_prediction_wall_s':None,
            'cost_missing_reason':'original interrupted parent saved fitted states but did not persist per-method timing; unknown, never zero-filled'})
        lineage.append({'method':tag,'original_state':fp(old/'GP_STATE.json'),'copied_state':fp(folder/'GP_STATE.json'),
            'original_prediction':fp(old/'TEST_PREDICTIONS_FROZEN.npz'),'selected_indices':indices,'new_prediction':fp(oldpred),'GP_refitted':False})
    rows=[r for r in all_rows if r['key'] in used]
    for row in rows:
        if not row['valid']:raise ValueError('invalid cache label')
        for field in ('result','receipt'):
            path=Path(row[field]['path'])
            if fp(path)['sha256']!=row[field]['sha256']:raise ValueError('actual query evidence changed')
        receipt=json.loads(Path(row['receipt']['path']).read_text())
        if receipt['failure'] is not None or receipt['code']['sha256']!=binding['source']['sha256']:raise ValueError('cache source/run receipt mismatch')
    validate_prediction_models(methods,rows,test,policy,out)
    (out/'physics').mkdir();write_json(out/'physics/REGISTRY.json',rows);write_json(out/'METHODS.json',methods)
    write_json(out/'QUERY_REQUESTS.json',[{'key':k,'role':'cache_import_model_membership_not_original_request_history','cache_hit':True} for k in sorted(used)])
    write_json(out/'CACHE_MIGRATION.json',{'parent_binding':binding,'parent_source':fp(snapshot),'new_source':fp(__file__),
        'same_physical_solver':solver_identity(),'lineage':lineage,'queries_retained':len(rows),'original_parent_completed':False,
        'test_labels_accessed':False,'new_PDE_queries':0,'per_method_timing':'not recoverable; null, not imputed'})
    def local_fp(p):
        p=Path(p);a=fp(p);a['path']=str(p.relative_to(out));return a
    seal={'phase':'training_complete_test_not_run','predictions':[local_fp(out/m['folder']/'TEST_PREDICTIONS_FROZEN.npz') for m in methods],
        'model_states':[local_fp(out/m['folder']/'GP_STATE.json') for m in methods],'policy':fp(policy_path),
        'test_design':fp(bound_file(policy_path,policy['test_design'])),'source':fp(__file__),'solver':solver_identity(),
        'methods':local_fp(out/'METHODS.json'),'registry':local_fp(out/'physics/REGISTRY.json'),'requests':local_fp(out/'QUERY_REQUESTS.json'),
        'terrain':fp(Path(args.inputs)/'copernicus_support.npz'),'numerics':{k:getattr(args,k) for k in ['order','backend','dry','cfl','maxdt']},
        'no_test_label_coordinates_in_registry':True,'migration':local_fp(out/'CACHE_MIGRATION.json')}
    if any(np.min(np.linalg.norm(np.asarray(r['x'])[None,:]-test,axis=1))<1e-12 for r in rows):raise ValueError('test/training collision')
    write_json(out/'PHASE1_SEAL.json',seal)


def run_formal_test(out,args):

    training=Path(args.training_dir);seal_path=training/'PHASE1_SEAL.json';seal=json.loads(seal_path.read_text());policy_path=Path(args.policy);policy=json.loads(policy_path.read_text())
    expected=len(policy['seeds'])*len(policy['methods'])
    if seal['phase']!='training_complete_test_not_run' or len(seal['predictions'])!=expected or len(seal['model_states'])!=expected:raise ValueError('complete all-model phase1 seal required')
    if fp(__file__)['sha256']!=seal['source']['sha256'] or fp(policy_path)['sha256']!=seal['policy']['sha256']:raise ValueError('test source/policy changed since training seal')
    if solver_identity()['sha256']!=seal['solver']['sha256']:raise ValueError('test solver identity differs')
    validate_execution(policy,args,policy_path)
    if fp(Path(args.inputs)/'copernicus_support.npz')['sha256']!=seal['terrain']['sha256']:raise ValueError('test terrain differs')
    for key,value in seal['numerics'].items():
        if getattr(args,key)!=value:raise ValueError('test numerical settings differ: '+key)
    for record in seal['predictions']+seal['model_states']+[seal['methods'],seal['registry'],seal['requests']]:
        p=training/record['path']
        if not p.is_file() or fp(p)['sha256']!=record['sha256']:raise ValueError('sealed training artifact changed: '+str(p))
    test_path=bound_file(policy_path,policy['test_design']);td=json.loads(test_path.read_text());test=np.r_[np.array(td['global']),np.array(td['development_predicted_boundary'])]
    methods=json.loads((training/'METHODS.json').read_text());old_records=json.loads((training/'physics/REGISTRY.json').read_text())
    validate_prediction_models(methods,old_records,test,policy,training)
    if any(np.min(np.linalg.norm(np.array(r['x'])[None,:]-test,axis=1))<1e-12 for r in old_records):raise ValueError('test truth already present in sealed training library')
    for p in training.rglob('TEST_KEYS.json'):raise ValueError('training directory already contains test labels')
    access=training/'TEST_ACCESS_STARTED.json'
    with access.open('x') as f:
        json.dump({'phase1_seal':fp(seal_path),'source':fp(__file__),'test_output':str(out),'purpose':'first independent test access; not training/tuning'},f,ensure_ascii=False,indent=2)
    write_json(out/'PHASE2_GATE.json',{'phase1_seal':fp(seal_path),'source':fp(__file__),'policy':fp(policy_path),'test_design':fp(test_path),'all_checks':'pass','test_queries_before_gate':0})
    args.horizon=policy['horizon_s'];registry=PhysicsRegistry(out/'physics',args,policy);registry.records={r['key']:r for r in old_records}
    test_high=registry.batch(test,'H','shared_independent_test_H24');test_low=registry.batch(test,'L','shared_independent_test_L24')
    write_json(out/'METHODS.json',methods);write_json(out/'TEST_KEYS.json',{'H':test_high,'L':test_low});write_json(out/'QUERY_REQUESTS.json',registry.requests)
    evaluate_formal_results(out,policy,methods,registry,test_high,test_low,training)


def evaluate_formal_results(out,policy,methods,registry,test_high,test_low,training_root):
    import csv
    rowsH=[registry.records[k] for k in test_high];rowsL=[registry.records[k] for k in test_low]
    truth=np.array([r.get('response_g') if r['valid'] else np.nan for r in rowsH],dtype=float)
    low=np.array([r.get('response_g') if r['valid'] else np.nan for r in rowsL],dtype=float)
    numeric_gray=np.isfinite(truth)&np.isfinite(low)&(abs(truth)<=abs(truth-low))
    rows=[];costs=[]
    for method in methods:
        folder=Path(training_root)/method['folder'];f=np.load(folder/'TEST_PREDICTIONS_FROZEN.npz');pred=f['mean_g']
        keys=list(dict.fromkeys(method['low_keys']+method['high_keys']));q=[registry.records[k] for k in keys]
        costs.append({'seed':method['seed'],'method':method['method'],'LF_queries':len(method['low_keys']),'HF_queries':len(method['high_keys']),
            'physical_child_wall_sum_s':sum(r['total_wall_s'] for r in q),'solver_wall_sum_s':sum(r['solver_wall_s'] for r in q),
            'compile_load_sum_s':sum(r['compile_and_load_s'] for r in q),'GP_fit_wall_s':method['fit_wall_s'],'acquisition_wall_s':method['acquisition_wall_s'],'auxiliary_prediction_wall_s':method['auxiliary_prediction_wall_s'],
            'logical_cold_start_total_s':sum(r['total_wall_s'] for r in q)+sum(method[k] for k in ['fit_wall_s','acquisition_wall_s','auxiliary_prediction_wall_s']) if all(method[k] is not None for k in ['fit_wall_s','acquisition_wall_s','auxiliary_prediction_wall_s']) else None})
        for label,index in [('global',np.arange(policy['test_global_count'])),('development_predicted_boundary',np.arange(policy['test_global_count'],len(truth)))]:
            valid=index[np.isfinite(truth[index])];y=truth[valid];p=pred[valid];actual=np.array([rowsH[i]['status']=='arrived' for i in valid]);positive=p>=0
            TP=int((actual&positive).sum());FN=int((actual&~positive).sum());FP=int((~actual&positive).sum());TN=int((~actual&~positive).sum())
            margin_true=.1*(2*np.exp(y)-1)-.1;margin_pred=np.maximum(0,.1*(2*np.exp(p)-1))-.1
            rows.append({'seed':method['seed'],'method':method['method'],'stratum':label,'n':len(index),'valid_n':len(valid),'unknown_n':len(index)-len(valid),
                'TP':TP,'FN':FN,'FP':FP,'TN':TN,'disagreement_rate':float((FN+FP)/len(valid)) if len(valid) else None,
                'arrival_miss_rate':FN/(TP+FN) if TP+FN else None,'RMSE_g':float(np.sqrt(np.mean((p-y)**2))) if len(valid) else None,
                'MAE_g':float(np.mean(abs(p-y))) if len(valid) else None,'RMSE_raw_margin_m':float(np.sqrt(np.mean((margin_true-margin_pred)**2))) if len(valid) else None,
                'numerically_unresolved_n':int(numeric_gray[index].sum()),'true_arrival_n':TP+FN,'time_proxy':'not_modelled; physical arrival brackets separately retained'})
    for filename,data in [('COMPARISON_METRICS.csv',rows),('METHOD_COSTS.csv',costs)]:
        with (out/filename).open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    write_json(out/'TEST_NUMERICAL_STATUS.json',{'gH':[float(x) if np.isfinite(x) else None for x in truth],'gL':[float(x) if np.isfinite(x) else None for x in low],
        'numerically_unresolved':numeric_gray.tolist(),'rule':policy['numerically_unresolved_rule'],'meaning':'empirical discrepancy, not a certified continuum error bound'})
    write_json(out/'EXPERIMENT_STATUS.json',{'generated':True,'methods':len(methods),'test_points':len(test_high),'test_query_unique_H':len(set(test_high)),'test_query_unique_L':len(set(test_low)),
        'physical_unique_queries':len(registry.records),'independent_review':'pending','test_used_for_tuning':False,'natural_validation':False})


def load_fixed_gp(state):

    fit=state['fit'];theta=np.asarray(fit['theta']);noise=fit['jitter_normalized_variance']
    if state['method']=='A_single_fill':
        m=SingleGP();m.x=np.asarray(state['xH']);m.mean=fit['mean'];m.scale=fit['scale'];m.y=(np.asarray(state['yH'])-m.mean)/m.scale
        K=matern52(m.x,m.x,np.exp(theta[0]),np.exp(theta[1:]))+noise*np.eye(len(m.x))
    else:
        m=ARGP();m.xl=np.asarray(state['xL']);m.xh=np.asarray(state['xH']);m.nL=len(m.xl);m.dim=m.xl.shape[1]
        m.ml=fit['mean_L'];m.mh=fit['mean_H'];m.scale=fit['output_scale'];m.y=np.r_[(np.asarray(state['yL'])-m.ml)/m.scale,(np.asarray(state['yH'])-m.mh)/m.scale]
        def unpack(t):
            k=m.dim;return np.exp(t[0]),np.exp(t[1:k+1]),np.exp(t[k+1]),np.exp(t[k+2:2*k+2]),t[-1]
        m.unpack=unpack;al,ll,ad,ld,rho=unpack(theta)
        A=matern52(m.xl,m.xl,al,ll);B=rho*matern52(m.xl,m.xh,al,ll)
        D=rho*rho*matern52(m.xh,m.xh,al,ll)+matern52(m.xh,m.xh,ad,ld)
        K=np.block([[A,B],[B.T,D]])+noise*np.eye(len(m.y))
    m.theta=theta;m.noise=noise;m.fit_record=fit;m.cf=cho_factor(K,lower=True);m.alpha=cho_solve(m.cf,m.y)
    return m


def cache_logical_id(point,level,args,policy):
    config=config_from_plane(point,policy);dx=policy['low_dx_m'] if level=='L' else policy['high_dx_m']
    body={'config':config,'dx':dx,'horizon':args.horizon,'dry':args.dry,'cfl':args.cfl,'backend':args.backend,
          'order':args.order,'maxdt':args.maxdt,'entry_width':args.entry_width,'solver':solver_identity()['sha256'],
          'assets':[(r['logical_id'],r['sha256']) for r in policy['input_assets']]}
    return hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest()


class FollowupRegistry(PhysicsRegistry):
    def __init__(self,out,args,policy):
        super().__init__(out,args,policy);self.aliases={};self.imported=[]
    def key(self,x,level):
        logical=cache_logical_id(x,level,self.args,self.policy)
        return self.aliases.get(logical,logical),config_from_plane(x,self.policy),self.policy['low_dx_m'] if level=='L' else self.policy['high_dx_m']
    def import_cache(self,path):
        if path is None:return
        parent=Path(path);rows=json.loads(parent.read_text());snapshot=parent.parent.parent/'SOURCE_SNAPSHOT.py'
        import importlib.util,sys
        spec=importlib.util.spec_from_file_location('followup_cache_source',snapshot);prior=importlib.util.module_from_spec(spec);sys.modules[spec.name]=prior;spec.loader.exec_module(prior)
        if prior.solver_identity()!=solver_identity():raise ValueError('cache physical solver changed')
        for row in rows:
            if not row['valid']:raise ValueError('invalid cached label')
            expected=config_from_plane(row['x'],self.policy)
            if row['config']!=expected or row['dx_m']!=(40 if row['level']=='L' else 30):raise ValueError('cached coordinate/config/fidelity mismatch')
            for field in ('receipt','result'):
                if fp(Path(row[field]['path']))['sha256']!=row[field]['sha256']:raise ValueError('cached evidence hash mismatch')
            receipt=json.loads(Path(row['receipt']['path']).read_text())
            if receipt['failure'] is not None or not receipt['source_unchanged_during_run'] or receipt['code']['sha256']!=fp(snapshot)['sha256']:raise ValueError('cache source or actual run mismatch')
            argv=receipt['argv']
            for flag,value in [('--horizon',self.args.horizon),('--dry',self.args.dry),('--cfl',self.args.cfl),('--order',self.args.order),('--maxdt',self.args.maxdt),('--entry-width',self.args.entry_width)]:
                if flag not in argv or float(argv[argv.index(flag)+1])!=value:raise ValueError('cache numerics mismatch: '+flag)
            if argv[argv.index('--backend')+1]!=self.args.backend or '--full' in argv:raise ValueError('cache execution scope mismatch')
            assets={Path(a['path']).name:a['sha256'] for a in receipt['inputs']}
            if any(assets.get(a['logical_id'])!=a['sha256'] for a in self.policy['input_assets']):raise ValueError('cache terrain/input assets changed')
            logical=cache_logical_id(row['x'],row['level'],self.args,self.policy)
            if logical in self.aliases and self.aliases[logical]!=row['key']:raise ValueError('conflicting logical cache entries')
            self.aliases[logical]=row['key'];self.records[row['key']]=row
            self.imported.append({'logical_id':logical,'original_key':row['key'],'original_receipt':row['receipt'],'original_result':row['result'],'origin_source':receipt['code'],'reused_without_PDE':True})
        write_json(self.out/'REGISTRY.json',list(self.records.values()));write_json(self.out/'CACHE_IMPORT.json',{'cache_registry':fp(parent),'source':fp(snapshot),'same_solver':solver_identity(),'entries':self.imported})


def import_seed137(out,args,registry,test):
    import shutil
    parent=Path(args.training_dir);seal=json.loads((parent/'PHASE1_SEAL.json').read_text());metas=json.loads((parent/'METHODS.json').read_text());lineage=[]
    if seal['solver']!=solver_identity():raise ValueError('seed137 model science differs')
    for meta in metas:
        if meta['seed']!=137:raise ValueError('only seed137 imported')
        relative=Path(meta['folder']);old=parent/relative;folder=out/relative;shutil.copytree(old,folder)
        for record in seal['model_states']:
            if Path(record['path'])==relative/'GP_STATE.json' and fp(old/'GP_STATE.json')['sha256']!=record['sha256']:raise ValueError('seed137 sealed state changed')
        state=json.loads((folder/'GP_STATE.json').read_text());model=load_fixed_gp(state);saved=np.load(old/'TEST_PREDICTIONS_FROZEN.npz');mean,std=model.predict(saved['x'])
        if not np.allclose(mean,saved['mean_g'],rtol=0,atol=1e-10) or not np.allclose(std,saved['std_g'],rtol=0,atol=1e-10):raise ValueError('fixed-parameter model restore differs')
        (folder/'TEST_PREDICTIONS_FROZEN.npz').rename(folder/'HISTORICAL24_PREDICTIONS.npz')
        aux=method_predictions(folder,model,test,meta['high_keys'],registry)
        meta={**meta,'fit_wall_s':None,'acquisition_wall_s':None,'auxiliary_prediction_wall_s':aux,'cost_missing_reason':'historical137 fit/acquisition missing; new548 prediction wall measured'}
        lineage.append({'method':meta['method'],'original_state':fp(old/'GP_STATE.json'),'new_state':fp(folder/'GP_STATE.json'),'GP_refitted':False,'historical_prediction_max_abs_difference':float(abs(mean-saved['mean_g']).max())})
        yield meta
    write_json(out/'SEED137_REUSE.json',{'parent_seal':fp(parent/'PHASE1_SEAL.json'),'entries':lineage,'new_PDE':0,'new_fit':0})


def run_followup_validation(out,args):

    training=Path(args.training_dir);policy_path=Path(args.policy);policy=json.loads(policy_path.read_text());seal=json.loads((training/'PHASE1_SEAL.json').read_text())
    validate_execution(policy,args,policy_path)
    if seal['source']['sha256']!=fp(__file__)['sha256'] or seal['policy']['sha256']!=fp(policy_path)['sha256']:raise ValueError('sealed source/policy changed')
    for rec in seal['predictions']+seal['model_states']+[seal['methods'],seal['registry'],seal['requests']]:
        if fp(training/rec['path'])['sha256']!=rec['sha256']:raise ValueError('sealed model/library changed')
    td=json.loads(bound_file(policy_path,policy['test_design']).read_text());test=np.r_[td['global'],td['development_predicted_boundary']]
    methods=json.loads((training/'METHODS.json').read_text());old=json.loads((training/'physics/REGISTRY.json').read_text());validate_prediction_models(methods,old,test,policy,training)
    with (training/'FOLLOWUP_TEST_ACCESS_STARTED.json').open('x') as f:json.dump({'source':fp(__file__),'seal':fp(training/'PHASE1_SEAL.json'),'labels_before_access':0},f)
    write_json(out/'PHASE2_GATE.json',{'source':fp(__file__),'seal':fp(training/'PHASE1_SEAL.json'),'all_declared_models_and532_predictions_immutable':True})
    registry=FollowupRegistry(out/'physics',args,policy)

    globalH=registry.batch(test[:16],'H','new_global16_H') if policy.get('run_new_global_truth',False) else []
    globalL=registry.batch(test[:16],'L','new_global16_L') if policy.get('run_new_global_truth',False) else []
    states={cut:{'evaluated':{},'sequence':[],'status':'initial'} for cut in (1,3)}
    def save_row(cut,j,key):
        r=registry.records[key];states[cut]['evaluated'][j]=key;states[cut]['sequence'].append({'j':j,'prediction_index':16+cut*129+j,'key':key,'status':r['status'],'h_m':r.get('target_max_h_m'),'valid':r['valid']})
    initial=[(cut,j) for cut in (1,3) for j in (0,64,128)]
    keys=registry.batch(np.array([test[16+c*129+j] for c,j in initial]),'H','critical_registered_initial6')
    for (cut,j),key in zip(initial,keys):save_row(cut,j,key)
    for cut,st in states.items():
        valid=all(registry.records[k]['valid'] for k in st['evaluated'].values());js=sorted(st['evaluated'])
        transitions=[(a,b) for a,b in zip(js[:-1],js[1:]) if (registry.records[st['evaluated'][a]]['status']=='arrived')!=(registry.records[st['evaluated'][b]]['status']=='arrived')] if valid else []
        st['initial_transitions']=transitions;st['status']='unknown_initial' if not valid else 'no_bracket' if not transitions else 'bracket_found'
        if transitions:st['lo'],st['hi']=transitions[0];st['left_positive']=registry.records[st['evaluated'][st['lo']]]['status']=='arrived'
    for iteration in range(4):
        planned=[(cut,(st['lo']+st['hi'])//2) for cut,st in states.items() if st['status']=='bracket_found']
        if not planned:break
        keys=registry.batch(np.array([test[16+c*129+j] for c,j in planned]),'H',f'critical_bisection_round{iteration+1}')
        for (cut,j),key in zip(planned,keys):
            save_row(cut,j,key);r=registry.records[key];st=states[cut]
            if not r['valid']:st['status']='unknown_bisection'
            elif (r['status']=='arrived')==st['left_positive']:st['lo']=j
            else:st['hi']=j
    searches=[];endpointH=[];endpointL=[];endpoint_indices=[]
    for cut,st in states.items():
        final=[st['lo'],st['hi']] if st['status']=='bracket_found' else None
        if final:
            st['status']='multiple_initial_transitions' if len(st['initial_transitions'])>1 else 'final_bracket'
            for j in final:endpointH.append(st['evaluated'][j]);endpoint_indices.append(16+cut*129+j)
        zv=(.125,.375,.625,.875)[cut]
        searches.append({'slice':cut,'zV':zv,'volume_m3':2e5*4**zv,'initial_transitions':st['initial_transitions'],'final_j':final,'status':st['status'],'HF_queries':len(st['sequence']),'sequence':st['sequence'],'final_mu_interval':[.03+.12*j/128 for j in final] if final else None,'discontinuity_warning':'labels bracket a discrete change; not a smooth h crossing'})
    if endpoint_indices:endpointL=registry.batch(test[endpoint_indices],'L','critical_final_endpoints_L') if policy.get('run_endpoint_low',False) else [None]*len(endpoint_indices)
    write_json(out/'CRITICAL_SEARCH.json',searches)
    write_json(out/'QUERY_REQUESTS.json',registry.requests);write_json(out/'TEST_KEYS.json',{'global_H':globalH,'global_L':globalL,'endpoint_H':endpointH,'endpoint_L':endpointL,'endpoint_prediction_indices':endpoint_indices})
    evaluate_followup(out,policy,methods,registry,training,globalH,globalL,endpointH,endpointL,endpoint_indices,searches)


def evaluate_followup(out,policy,methods,registry,training,gH,gL,eH,eL,indices,searches):
    import csv
    rows=[];boundary=[]
    for meta in methods:
        f=np.load(training/meta['folder']/'TEST_PREDICTIONS_FROZEN.npz')
        for label,H,L,ix in [('new_global',gH,gL,list(range(16))),('real_critical_endpoints',eH,eL,indices)]:
            good=[(h,l,i) for h,l,i in zip(H,L,ix) if registry.records[h]['valid']];y=np.array([registry.records[h]['response_g'] for h,l,i in good]);p=np.array([f['mean_g'][i] for h,l,i in good]);actual=np.array([registry.records[h]['status']=='arrived' for h,l,i in good],dtype=bool);pred=p>=0
            TP=int((actual&pred).sum());FN=int((actual&~pred).sum());FP=int((~actual&pred).sum());TN=int((~actual&~pred).sum())
            rows.append({'seed':meta['seed'],'method':meta['method'],'stratum':label,'n':len(H),'valid_n':len(good),'unknown_n':len(H)-len(good),'TP':TP,'FN':FN,'FP':FP,'TN':TN,'disagreement_rate':float(np.mean(actual!=pred)) if len(good) else None,'arrival_miss_rate':FN/(TP+FN) if TP+FN else None,'RMSE_g':float(np.sqrt(np.mean((p-y)**2))) if len(good) else None,'RMSE_raw_margin_m':float(np.sqrt(np.mean((.2*np.exp(y)-.1-np.maximum(0,.2*np.exp(p)-.1))**2))) if len(good) else None})
        for search in searches:
            cut=search['slice'];means=f['mean_g'][16+cut*129:16+(cut+1)*129];js=np.flatnonzero((means[:-1]>=0)!=(means[1:]>=0));roots=[]
            for j in js:roots.append(.03+.12*(j+abs(means[j])/(abs(means[j])+abs(means[j+1])))/128)
            interval=search['final_mu_interval'];distance=max(interval[0]-roots[0],0,roots[0]-interval[1]) if roots and interval else None
            boundary.append({'seed':meta['seed'],'method':meta['method'],'slice':cut,'true_bracket':interval,'predicted_zero_crossings':roots,'mu_boundary_distance_to_true_bracket':distance,'prediction_has_no_crossing':not roots,'search_status':search['status']})
    with (out/'COMPARISON_METRICS.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write_json(out/'BOUNDARY_LOCATION.json',boundary)
    write_json(out/'VALIDATION_COSTS.json',{'global_H':len(gH),'global_L':len(gL),'all_critical_search_H':sum(s['HF_queries'] for s in searches),'critical_endpoint_L':sum(k is not None for k in eL),'actual_unique_queries':len(registry.records),'child_wall_sum_s':sum(r['total_wall_s'] for r in registry.records.values()),'search_selection_cost_included':True})
    write_json(out/'EXPERIMENT_STATUS.json',{'methods':len(methods),'seeds':policy['seeds'],'new_global_coordinates':16,'new_global_points_actually_queried':len(gH),'critical_potential_coordinates':516,'actual_critical_H':sum(s['HF_queries'] for s in searches),'actual_endpoint_H':len(eH),'actual_endpoint_pairs':sum(k is not None for k in eL),'historical24_used_as_new_test':False,'test_returned_to_training':False,'independent_review':'pending'})


def seal_followup_subset(out,args):

    import shutil
    parent=Path(args.training_dir);policy_path=Path(args.policy);policy=json.loads(policy_path.read_text());validate_execution(policy,args,policy_path)
    if policy['seeds']!=[137]:raise ValueError('deadline137 subset only')
    methods=[m for m in json.loads((parent/'METHODS.json').read_text()) if m['seed']==137]
    if len(methods)!=4:raise ValueError('all four complete137 methods required')
    td=json.loads(bound_file(policy_path,policy['test_design']).read_text());test=np.r_[td['global'],td['development_predicted_boundary']];keys=set();lineage=[]
    for meta in methods:
        old=parent/meta['folder'];folder=out/meta['folder'];shutil.copytree(old,folder);keys.update(meta['low_keys']+meta['high_keys'])
        lineage.append({'method':meta['method'],'origin_state':fp(old/'GP_STATE.json'),'state':fp(folder/'GP_STATE.json'),'origin_prediction':fp(old/'TEST_PREDICTIONS_FROZEN.npz'),'prediction':fp(folder/'TEST_PREDICTIONS_FROZEN.npz'),'refit':False})
    rows=[r for r in json.loads((parent/'physics/REGISTRY.json').read_text()) if r['key'] in keys]
    validate_prediction_models(methods,rows,test,policy,out);(out/'physics').mkdir();write_json(out/'physics/REGISTRY.json',rows);write_json(out/'METHODS.json',methods);write_json(out/'QUERY_REQUESTS.json',[])
    write_json(out/'SUBSET_MIGRATION.json',{'parent_source':fp(parent/'SOURCE_SNAPSHOT.py'),'source':fp(__file__),'model_seeds':[137],'excluded_partial_seeds':[271,503],'lineage':lineage,'new_PDE':0,'refit':0,'reason':'18:40 deadline, scope selected before critical truth'})
    def local(p):
        p=Path(p);r=fp(p);r['path']=str(p.relative_to(out));return r
    write_json(out/'PHASE1_SEAL.json',{'phase':'training_complete_test_not_run','source':fp(__file__),'policy':fp(policy_path),'solver':solver_identity(),
        'predictions':[local(out/m['folder']/'TEST_PREDICTIONS_FROZEN.npz') for m in methods],'model_states':[local(out/m['folder']/'GP_STATE.json') for m in methods],
        'methods':local(out/'METHODS.json'),'registry':local(out/'physics/REGISTRY.json'),'requests':local(out/'QUERY_REQUESTS.json'),'migration':local(out/'SUBSET_MIGRATION.json')})


def decision_request(args):

    if args.full or args.dx!=30:raise ValueError('decision publication supports local30m only')
    config={'volume_m3':float(args.volume),'duration_s':float(args.duration),'mu':float(args.mu),'xi_m_s2':float(args.xi),
            'speed_m_s':float(args.speed),'entry_width_m':float(args.entry_width),'shape':args.shape,'terrain_name':args.terrain_name}
    if not all(np.isfinite(config[k]) for k in ['volume_m3','duration_s','mu','xi_m_s2','speed_m_s','entry_width_m']):raise ValueError('finite request required')
    body={'config':config,'numerics':{**{k:float(getattr(args,k)) for k in ['dx','horizon','cfl','dry','maxdt']},'order':int(args.order),'backend':str(args.backend)},
          'inputs':[{'logical_id':n,'sha256':fp(Path(args.inputs)/n)['sha256']} for n in ['copernicus_support.npz','aster_support.npz','TERRAIN.json','mapped_path_utm.csv']],
          'science':solver_identity(),'target':{'axis':'x','coordinate_m':348720,'tangent_interval_m':[3134880,3135600],
          'criterion':'instantaneous max vertical event depth >= threshold on fixed interpolated line','threshold_m':.1,'horizon_s':float(args.horizon),'full_reach':False}}
    return {'request_id':hashlib.sha256(json.dumps(body,sort_keys=True).encode()).hexdigest(),**body}


def physical_evidence_check(request,result_path,receipt_path,source_certificate):

    d=json.loads(Path(result_path).read_text());receipt=json.loads(Path(receipt_path).read_text())
    if receipt['failure'] is not None or not receipt['source_unchanged_during_run']:raise ValueError('failed or changed-source physical run')
    if receipt['code']['sha256']!=source_certificate['source']['sha256'] or source_certificate['science']!=request['science']:raise ValueError('physical scientific identity mismatch')
    if not any(Path(r['path']).name==Path(result_path).name and r['sha256']==fp(Path(result_path))['sha256'] for r in receipt['outputs']):raise ValueError('result not bound by actual run receipt')
    inputs={Path(i['path']).name:i['sha256'] for i in receipt['inputs']}
    if any(inputs.get(i['logical_id'])!=i['sha256'] for i in request['inputs']):raise ValueError('physical input asset mismatch')
    argv=receipt['argv'];c=request['config'];n=request['numerics']
    numeric={'--volume':c['volume_m3'],'--duration':c['duration_s'],'--mu':c['mu'],'--xi':c['xi_m_s2'],'--speed':c['speed_m_s'],
             '--entry-width':c['entry_width_m'],'--dx':n['dx'],'--horizon':n['horizon'],'--cfl':n['cfl'],'--dry':n['dry'],'--order':n['order'],'--maxdt':n['maxdt']}
    for flag,value in numeric.items():
        if flag not in argv or float(argv[argv.index(flag)+1])!=value:raise ValueError('physical config/numerics mismatch: '+flag)
    for flag,value in [('--backend',n['backend']),('--shape',c['shape']),('--terrain-name',c['terrain_name'])]:
        if flag not in argv or argv[argv.index(flag)+1]!=value:raise ValueError('physical mode mismatch: '+flag)
    if '--full' in argv or d['full_reach']:raise ValueError('physical domain scope mismatch')
    support=d['arrival_support'];target=request['target']
    if support['axis']!=target['axis'] or support['coordinate_m']!=target['coordinate_m'] or support['tangent_interval_m']!=target['tangent_interval_m']:raise ValueError('physical target support mismatch')
    if not d['whole_window_valid'] or not np.isfinite(d['valid_until_s']) or d['valid_until_s']<n['horizon']:raise ValueError('domain/time window invalid; unknown')
    for field in ['mass_error_rel','momentum_error_rel']:
        if not np.isfinite(d[field]) or abs(d[field])>1e-8:raise ValueError('physical budget check failed: '+field)
    if not np.isfinite(d['target_max_h']) or d['target_max_h']<0 or not np.isfinite(d['min_h']) or d['min_h']<0:raise ValueError('invalid physical thickness')
    arrival=d['threshold_arrivals']['0.1'];status=arrival['status'];bracket=arrival['time_bracket_s']
    if status=='arrived':
        if d['target_max_h']<.1 or bracket is None or not 0<=bracket[0]<=bracket[1]<=n['horizon']:raise ValueError('arrival criterion/bracket inconsistent')
    elif status=='valid_no_arrival_by_horizon':
        if d['target_max_h']>=.1 or bracket is not None:raise ValueError('no-arrival criterion inconsistent')
    else:raise ValueError('no valid physical class; unknown')
    return {'status':status,'time_bracket_s':bracket,'target_max_h_m':d['target_max_h'],'result':fp(Path(result_path)),'receipt':fp(Path(receipt_path))}


def bind_decision_cache(out,args):

    import importlib.util,sys,shutil
    request=decision_request(args);spec=importlib.util.spec_from_file_location('guard_evidence_source',args.evidence_source);prior=importlib.util.module_from_spec(spec);sys.modules[spec.name]=prior;spec.loader.exec_module(prior)
    certificate={'source':fp(args.evidence_source),'science':prior.solver_identity(),'binding_source':fp(__file__)}
    checked=physical_evidence_check(request,args.evidence_result,args.evidence_receipt,certificate)
    evidence=out/'evidence';evidence.mkdir();shutil.copy2(args.evidence_result,evidence/'pilot_30.json');shutil.copy2(args.evidence_receipt,evidence/'RUN_RECEIPT.json')
    entry={'schema':'strict-physical-verification-cache-v1','request':request,'source_certificate':certificate,
           'result':{'package_path':'evidence/pilot_30.json','sha256':fp(evidence/'pilot_30.json')['sha256']},
           'receipt':{'package_path':'evidence/RUN_RECEIPT.json','sha256':fp(evidence/'RUN_RECEIPT.json')['sha256']},'bound_physical_status':checked['status']}
    write_json(out/'CACHE_ENTRY.json',entry)


def load_decision_cache(path,request):
    path=Path(path);entry=json.loads(path.read_text())
    if entry['schema']!='strict-physical-verification-cache-v1' or entry['request']!=request:raise ValueError('cache request identity mismatch')
    resolved=[]
    for field in ('result','receipt'):
        relative=Path(entry[field]['package_path']);p=path.parent/relative
        if relative.is_absolute() or not p.resolve().is_relative_to(path.parent.resolve()):raise ValueError('cache package path unsafe')
        if not p.is_file() or fp(p)['sha256']!=entry[field]['sha256']:raise ValueError('cache file absent or changed: '+field)
        resolved.append(p)
    return physical_evidence_check(request,*resolved,entry['source_certificate'])


def run_decision(out,args):
    import datetime
    started=time.perf_counter();request=decision_request(args);estimate={'status':'not_available','authoritative':False};reason=None;verified=None;kind='none';fresh_cost=None
    if args.gp_state is not None:
        state=json.loads(args.gp_state.read_text());c=request['config']
        in_scope=(2e5<=c['volume_m3']<=8e5 and .03<=c['mu']<=.15 and c['duration_s']==960 and c['xi_m_s2']==500 and c['speed_m_s']==30 and c['entry_width_m']==30 and c['terrain_name']=='copernicus' and c['shape']=='triangle' and args.horizon==1800 and args.cfl==.35 and args.dry==1e-4 and args.order==2 and args.maxdt==2)
        if in_scope:
            model=load_fixed_gp(state);x=np.array([[np.log(c['volume_m3']/2e5)/np.log(4),(c['mu']-.03)/.12]]);m,s=model.predict(x)
            estimate={'status':'arrived' if m[0]>=0 else 'valid_no_arrival_by_horizon','mean_g':float(m[0]),'std_g':float(s[0]),'model_state':fp(args.gp_state),'authoritative':False}
        else:estimate={'status':'outside_frozen_GP_scope','authoritative':False,'model_state':fp(args.gp_state)}
    write_json(out/'REQUEST_AND_GP_FROZEN.json',{'frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':fp(__file__),'request':request,'GP_estimate':estimate,'physical_truth_accessed_yet':False})
    if args.cache_entry is not None:
        try:verified=load_decision_cache(args.cache_entry,request);kind='exact_fingerprint_cache'
        except Exception as exc:reason=type(exc).__name__+": "+str(exc)
    if verified is None and not args.cache_only:
        try:
            row=execute_physics(out/'physical_verification',request['config'],args,30.)
            write_json(out/'PHYSICAL_CALL.json',row);fresh_cost={'child_wall_s':row['total_wall_s'],'solver_wall_s':row.get('solver_wall_s'),'compile_and_load_s':row.get('compile_and_load_s')}
            if row.get('valid'):
                try:verified=physical_evidence_check(request,row['result']['path'],row['receipt']['path'],{'source':fp(__file__),'science':solver_identity()});kind='fresh_30m_PDE'
                except (ValueError,KeyError,OSError,TypeError) as exc:reason=str(exc)
            else:reason='physical query failed or invalid domain/window; no class published'
        except Exception as exc:
            reason="physical verification exception: "+type(exc).__name__+": "+str(exc)
            write_json(out/"VERIFICATION_EXCEPTION.json",{"error_type":type(exc).__name__,"message":str(exc),"GP_fallback_forbidden":True})
    final={'request':request,'GP_estimate':estimate,'final_status':verified['status'] if verified else 'unknown','needs_verification':verified is None,
           'arrival_time_bracket_s':verified['time_bracket_s'] if verified else None,'verification_source':kind,'physical_evidence':verified,
           'unverified_reason':reason if not verified else None,'fresh_physical_cost':fresh_cost,'decision_wall_s':time.perf_counter()-started,
           'published_from_GP':False,'scope':'D_ONLY conditional physical-model verification; not natural-event zero-error guarantee'}
    write_json(out/'DECISION.json',final)


def sol():
    import argparse,sys,platform,datetime,importlib.metadata
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['bench','terrain','pilot','smooth','morris','freeze','near','perturb','train','test','global','inverse','baseline','join','seal-cache','validate-new','seal-subset','decide','bind-cache']);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--inputs',type=Path);ap.add_argument('--dx',type=float,default=30);ap.add_argument('--dry',type=float,default=1e-4);ap.add_argument('--order',type=int,choices=[1,2],default=2);ap.add_argument('--backend',choices=['numpy','c'],default='numpy');ap.add_argument('--workers',type=int,choices=[1,2,3],default=2);ap.add_argument('--morris-dir',type=Path);ap.add_argument('--policy',type=Path);ap.add_argument('--training-dir',type=Path);ap.add_argument('--previous-design',type=Path);ap.add_argument('--full-inputs',type=Path);ap.add_argument('--q3-doc',type=Path);ap.add_argument('--q4-doc',type=Path);ap.add_argument('--entity-root',type=Path);ap.add_argument('--allow-synthetic',action='store_true')
    ap.add_argument('--gp-state',type=Path);ap.add_argument('--cache-entry',type=Path);ap.add_argument('--cache-only',action='store_true');ap.add_argument('--evidence-result',type=Path);ap.add_argument('--evidence-receipt',type=Path);ap.add_argument('--evidence-source',type=Path);ap.add_argument('--cache-registry',type=Path);ap.add_argument('--entry-width',type=float,default=30.);ap.add_argument('--full',action='store_true');ap.add_argument('--volume',type=float,default=8e5);ap.add_argument('--duration',type=float,default=960);ap.add_argument('--speed',type=float,default=30);ap.add_argument('--mu',type=float,default=.05);ap.add_argument('--xi',type=float,default=500);ap.add_argument('--horizon',type=float,default=1800);ap.add_argument('--cfl',type=float,default=.35);ap.add_argument('--maxdt',type=float,default=2);ap.add_argument('--terrain-name',default='copernicus',choices=['copernicus','aster']);ap.add_argument('--shape',default='triangle',choices=['triangle','double']);args=ap.parse_args()
    global DRY,ORDER,FAST
    DRY=args.dry;ORDER=args.order
    args.out.mkdir(parents=True,exist_ok=False);start=datetime.datetime.now(datetime.timezone.utc).isoformat();code_before=fp(__file__)
    (args.out/'SOURCE_SNAPSHOT.py').write_bytes(Path(__file__).read_bytes())
    root=Path(__file__).resolve().parents[3] if 'work/q4/' in str(Path(__file__).resolve()) else Path.cwd()

    compile_metadata=None
    inputs=[];error=None
    try:
        if args.inputs:inputs=[fp(p) for p in sorted(args.inputs.iterdir()) if p.is_file()]
        if args.backend=='c' and args.mode not in ('morris','freeze','near','perturb','train','test','global','inverse','baseline','join','seal-cache','validate-new','seal-subset','decide','bind-cache','terrain'):
            tic=time.perf_counter();FAST=FastFlux(args.out/'cbuild');FAST.dry=DRY
            compile_metadata={'C_source_sha256':FAST.sha,'library_sha256':FAST.library_sha256,'compile_and_load_s':time.perf_counter()-tic,'actual_compile_command':FAST.command,'fast_math':False}
        if args.mode=='bench':benchmarks(args.out)
        elif args.mode=='smooth':smooth_benchmark(args.out)
        elif args.mode=='morris':run_morris(args.out,args)
        elif args.mode=='freeze':freeze_formal_policy(args.out,args)
        elif args.mode=='near':run_near_boundary(args.out,args)
        elif args.mode=='perturb':run_perturbations(args.out,args)
        elif args.mode=='train':run_formal_training(args.out,args)
        elif args.mode=='test':run_formal_test(args.out,args)
        elif args.mode=='decide':run_decision(args.out,args)
        elif args.mode=='bind-cache':bind_decision_cache(args.out,args)
        elif args.mode=='seal-subset':seal_followup_subset(args.out,args)
        elif args.mode=='validate-new':run_followup_validation(args.out,args)
        elif args.mode=='seal-cache':seal_cached_seed(args.out,args)
        elif args.mode=='global':run_global_cases(args.out,args)
        elif args.mode=='inverse':run_synthetic_inverse(args.out,args)
        elif args.mode=='baseline':write_json(args.out/'AXIS_BASELINE.json',along_axis_baseline(args.inputs,args.full,args.mu,args.xi,horizon=args.horizon))
        elif args.mode=='join':write_json(args.out/'Q3_Q4_JOIN.json',strict_q3_q4_join(json.loads(args.q3_doc.read_text()),json.loads(args.q4_doc.read_text()),args.entity_root,args.allow_synthetic))
        elif args.mode=='terrain':terrain(args.out,root,args.full)
        else:pilot(args.out,args.inputs,args.dx,args)
    except Exception as exc:
        error={'type':type(exc).__name__,'message':str(exc)};raise
    finally:
        write_json(args.out/'RUN_RECEIPT.json',{'started_at':start,'finished_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'argv':sys.argv,'code':code_before,'code_after':fp(__file__),'acceleration':compile_metadata,'source_unchanged_during_run':code_before==fp(__file__),'inputs':inputs,'environment':{'python':sys.version,'platform':platform.platform(),'versions':runtime_versions()},'failure':error,'outputs':[fp(p) for p in sorted(args.out.iterdir()) if p.is_file() and p.name!='RUN_RECEIPT.json'],'scientific_approval':'pending'})


if __name__=='__main__':sol()


