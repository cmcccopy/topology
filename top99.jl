"""拓扑优化99代码Julia版"""
using LinearAlgebra
using SparseArrays
using Plots

"""建立单元刚度矩阵\n
    平面4节点单元,每节点2自由度
"""
function KE()
    μ=0.3   #泊松比
    a=1.;b=1.;
    E=1.#弹性模量
    t=1.#单元厚度
    H=1/(1-μ^2);   r=(1-μ)/2;  s=(1-3μ)/8;  
    m=(1+μ)/8;     α=a/(3b);      β=b/(3a)
    ke=H*E*t*[β+r*α      0 0 0 0 0 0 0;
              m          α+β*r    0 0 0 0 0 0;
             -β+0.5*r*α  s          β+r*α 0 0 0 0 0;
             -s          α/2-r*β    -m           α+β*r 0 0 0 0;
             (-β-r*α)/2  -m         β/2-r*α      s             β+r*α 0 0 0;
             -m          -α/2-r*β/2  -s           -α+(r*β)/2    m           α+β*r 0 0;
             β/2-r*α     -s         (-β-r*α)/2   m             -β+r*α/2    s          β+r*α 0;
             s           -α+(r*β)/2  m           (-α-r*β)/2    -s          α/2-r*β    -m     α+β*r]
   ke=Symmetric(ke,:L)#以：L形式形成对称矩阵
    return ke
end

"""有限元求解程序"""
function fem(nelx::Int,nely::Int,x::Matrix{Float64},penal::Number)
  
    ke=KE()
    K=spzeros(2*(nelx+1)*(nely+1),2*(nelx+1)*(nely+1))#创建整体单元矩阵，以稀疏矩阵形式
    F = zeros(2*(nely+1)*(nelx+1)); 
    U = zeros(2*(nely+1)*(nelx+1));

    #"""组装整体刚度矩阵"""
    for elx = 1:nelx
        for ely = 1:nely
          n1 = (nely+1)*(elx-1)+ely; 
          n2 = (nely+1)* elx   +ely;
          edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
          K[edof,edof] = K[edof,edof] + x[ely,elx]^penal*ke
        end
      end
   
    F[2]=-1#向下的力
    fixdofs=union(1:2:2*(nely+1),2*(nelx+1)*(nely+1))#固定的自由度标号
    for i=fixdofs #对固定自由度的整体刚度矩阵处理
      K[i,i]=K[i,i]*10e8
      F[i]=0
    end
    U=K\F
    
  return U
end

"""灵敏度过滤程序

"""
function check(nelx::Int64, nely::Int64, rmin::Number, x::Matrix{Float64}, dc::Matrix{Float64})
  dcn=zeros(nely,nelx)#灵敏度初始化
  for i=1:nelx
    for j=1:nely 
      sum=0.
      for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)        #floor为向下取整函数，遍历半径rmin内的单元
        for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)      #得到ij目标单元xy两方向的单元数并遍历
          fac = rmin-sqrt((i-k)^2+(j-l)^2);
          sum = sum+max(0,fac);
          dcn[j,i] = dcn[j,i] + max(0,fac)*x[l,k]*dc[l,k];
        end
      end
      dcn[j,i]=dcn[j,i]/(x[j,i]*sum)
    end
  end
  return dcn
end

"""OC优化法求变量x"""
function OC(nelx,nely,x,volfrac,dc)  
  l1 = 0; l2 = 100000; move = 0.2;
  xnew = zeros(nely,nelx);    
  while l2-l1 > 1e-4 #用二分法求λ₁
      lmid = 0.5*(l2+l1);
      Cᵢ = sqrt.(-dc/lmid) ; 
      xᵏCᵢ = x.*Cᵢ ;   
      xnew= max.(fill(0.001,size(x)),max.(x.-move,min.(fill(1.,size(x)),min.(x.+move,xᵏCᵢ))));
      if sum(sum(xnew)) - volfrac*nelx*nely > 0;
          l1 = lmid;
      else
          l2 = lmid; 
      end
  end
  return xnew
end 
 


"""主程序"""
function top(nelx::Int64,nely::Int64,volfrac::Number,penal::Number,rmin::Number)
  x=fill(volfrac,(nely,nelx))#初始化密度x
    loop=0
    chang=1.
    dc=zeros(nely,nelx)
    C_Values=[];            #储存柔度C的值变化情况
    
    while chang>0.01
      loop+=1
      xold=x
      U=fem(nelx,nely,x,penal)
      ke=KE()
      C=0.
                #求解柔度C和灵敏度dc
        for ely=1:nely
          for elx=1:nelx 
            n1 = (nely+1)*(elx-1)+ely; 
            n2 = (nely+1)* elx+ely;
            edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2]
            Ue=U[edof];
            
            C=C+x[ely,elx]^penal*Ue'*ke*Ue

            dc[ely,elx]=-penal*x[ely,elx]^(penal-1)*Ue'*ke*Ue
          end
        end

        push!(C_Values,C)#将C存储到C_Values

        dc=check(nelx,nely,rmin,x,dc)
        
        x=OC(nelx,nely,x,volfrac,dc)

       chang=maximum(abs.(x-xold))#循环判定依据


    end
    return x,C_Values,loop
end






nelx=50;nely=33; #x方向，y方向单元个数
penal=3.         #x的惩罚因子
rmin=2           #网格灵敏度过滤半径
volfrac=0.5      #体积分数
x,c,l=top(nelx,nely,volfrac,penal,rmin)
plot(1:l,c,title = "sensity evolution test", label = "sensity c") 
xlabel!("iteration")
heatmap(x)  #绘制矩阵的热图



