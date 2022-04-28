using Flux,Statistics,Plots,MLDataUtils,DelimitedFiles
using Flux.Data: DataLoader
#provided there is an iscaler function the iscalerbatch and inv_iscaler function will respectively scale and unscale a matrix dataset
function iscalerbatch(x)
    nfeat,sl = size(x)
    out=Array{Float32,2}(undef,nfeat,sl)
    for i=1:sl
        out[:,i]=iscaler(x[:,i])
    end
    return out
end
function inv_iscalerbatch(x)
    nfeat,sl = size(x)
    out=Array{Float32,2}(undef,nfeat,sl)
    for i=1:sl
        out[:,i]=inv_iscaler(x[:,i])
    end
    return out
end
#Batchseqmat and unbatchseqmat will turn a n dimensional vector of nfeat x seqlen data into a seqlen vector of nfeat x n and unbatchseq will perform the transform in reverse
function batchseqmat(X)
    #Here X is a vector of matrices
    nbatch = size(X)[1]
    nfeat,sl= size(X[1])
    out = [Array{Float32}(undef,nfeat,nbatch) for i=1:sl]
    for i=1:sl
        for j=1:nbatch
            out[i][:,j]=X[j][:,i]
        end
    end
    return out
end
function unbatchseqmat(X)
    sl = size(X,1)
    nfeat,nbatch= size(X[1])
    out = [Array{Float32}(undef,nfeat,sl) for i=1:nbatch]
    for i=1:nbatch
        for j=1:sl
            out[i][:,j]=X[j][:,i]
        end
    end
    return out
end
#Calculates mean percent error on apriori test set
function mean_error_func(x,y,NN)
    yNN = NN(x);
    yNNcomp = inv_oscaler(yNN)[:]
    ycomp = inv_oscaler(y)[:]
    return mean(100*abs.(ycomp.-yNNcomp)./ycomp)
end
#calculates mean squared error using dynamic models
function cust_mse(x,y,NN)
    Flux.reset!(NN)
    return mean((NN(x)[:]-y[:]).^2)
end
#calculates root mean squared errors using dynamic models
function RMSE(x, y,NN)
     sqrt(cust_mse(x,y,NN))
end
#scaler and inverse scaler for temporal dimension
tscale(t)=Float32(t/(300/3600));
inv_tscale(t_sc)=Float32(t_sc*(300/3600));
#evalprofile will use a dynamic model m to construct a 20 node discrete accumulation profile
#A recurrent model rm will require you to perform Flux.reset!(rm) prior to using this function
function evalprofile(s,m)
    s=reshape(s,(1,3))
    x=[s 0f0 0f0]
    x = iscaler(Float32.(reshape(x,(nfeat,1))))
    y = Array{Float32,2}(undef,sl+1,1)
    y[1]=0
    t_i=0
    y_i=0
    for i=2:sl+1
        y_i=m(x)[1]
        y[i]=y_i
        x[4]+=.05;
        x[5]=y_i;
        x, y_i
    end
    return inv_oscaler(y)[2:end]
end
#Pe evaluates the percent error between a true dataset yt and a predicted dataset yp
pe(yt,yp)=100*abs.(yt.-yp)./yt
#EVALMPE builds dynamic profiles using only biological parameters and uses the constructed profiles to evaluate mean percent error
#This is equivalent to taking mean percent error using baseline multidimensional regression model
#evalmpe takes in a matrix S containing random samples of the biological parameters lp,k,rs the corresponding matrix YS containing 20 node accumulation profiles and a dynamic model NN
#eval mpe will then output a mean percent error between the model and the dataset by constructing the profile for each sample S[:,i]
function evalmpe(S,YS,nn)
    Flux.reset!(nn)
    n=size(S,2)
    out=Array{Float32,2}(undef,20,n)
    for i=1:n
        ypred=evalprofile(S[:,i],nn)
        out[:,i]=pe(YS[:,i],ypred)
    end
    return mean(out)
end
