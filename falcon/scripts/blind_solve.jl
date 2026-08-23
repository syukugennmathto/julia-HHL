#!/usr/bin/env julia
#
# blind_solve.jl -- recover the key from first-two events WITHOUT knowing which
# of the two calls straddled (the blind which-call disjunction of section 6.5).
#
#     julia --project=falcon falcon/scripts/blind_solve.jl
#
# Each event i gives a public syndrome c_i and the single bit "a first-two
# centre was an integer".  The two candidate rows are a_i = form_row(c_i, n/2-1)
# and b_i = form_row(c_i, n-1), and one satisfies row.f = 0 (mod q).  A blind
# attacker does not know which -- naive iterative solving (EM) does not resolve
# it, and no public statistic of the two signatures separates the two calls
# (both desync at the very start of the descent).
#
# The structure that cracks it: b_i = x^{n/2} a_i, so the disjunction is
#     (a_i . f)(a_i . f') = 0,   f' = -x^{n/2} f,
# a LINEAR measurement on the rank-2 symmetric matrix S = sym(f f'^T):
#     a_i^T S a_i = 0.
# With M ~ n(n+1)/2 events the measurements pin S (kernel dim 1); its column
# space is exactly span{f, f'}, from which f is recovered as the short vector.
# So the blind channel is solvable -- at O(n^2) events, an n-fold overhead over
# the O(n) events the labelled channel needs (section 6.5).
#
# This runs on synthetic events with real q (the algebra is identical at n=512;
# only the O(n^2) linear solve is why we demonstrate at n <= 32).
using Printf, Random
const q = 12289
sh(v,m)  = (n=length(v); Int[ (i=j-m; i>=0 ? v[i+1] : -v[i+n+1]) for j in 0:n-1 ])
frow(c,k,n) = Int[ (i=k-j; mod(i>=0 ? c[i+1] : -c[i+n+1], q)) for j in 0:n-1 ]
inv_(a)  = invmod(mod(a,q),q)
ctr(v)   = Int[(y=mod(c,q); y>q÷2 ? y-q : y) for c in v]
function rref(M0)
    M=[mod(x,q) for x in M0]; rows,cols=size(M); piv=Int[]; r=1
    for col in 1:cols
        pr=findfirst(i->M[i,col]%q!=0, r:rows); pr===nothing&&continue; pr+=r-1
        M[r,:],M[pr,:]=M[pr,:],M[r,:]; iv=inv_(M[r,col]); M[r,:]=mod.(M[r,:].*iv,q)
        for i in 1:rows; i==r&&continue; fc=M[i,col]; fc==0&&continue; M[i,:]=mod.(M[i,:].-fc.*M[r,:],q); end
        push!(piv,col); r+=1; r>rows&&break
    end
    M,piv
end
function kernelbasis(A,ncol)
    R,piv=rref(A); free=setdiff(1:ncol,piv); B=Vector{Int}[]
    for fc in free; v=zeros(Int,ncol); v[fc]=1; for (ri,c) in enumerate(piv); v[c]=mod(-R[ri,fc],q); end; push!(B,v); end
    B
end
symmeas(a,n)=(D=n*(n+1)÷2; out=zeros(Int,D); t=1; for i in 1:n,j in i:n; out[t]= i==j ? mod(a[i]*a[j],q) : mod(2*a[i]*a[j],q); t+=1; end; out)
tomat(vec,n)=(S=zeros(Int,n,n); t=1; for i in 1:n,j in i:n; S[i,j]=vec[t]; S[j,i]=vec[t]; t+=1; end; S)
# two independent columns of S = a basis of colspace(S)=span{f,f'}
function colbasis(S,n)
    b=Vector{Int}[]
    for j in 1:n
        c=S[:,j]; any(!=(0),c) || continue
        if isempty(b); push!(b,c)
        else
            # independent from b[1] mod q?
            R,p=rref(hcat(b[1],c)'|>Matrix); if length(p)==2; push!(b,c); break; end
        end
    end
    b
end
# short vector in span{g1,g2}: enumerate q+1 projective directions; anchor scale
# on the first nonzero coord being small (|f_i| is tiny)
function shortvec(g1,g2,n; bnd=8)
    best=nothing; bn=typemax(Int)
    dirs = vcat([mod.(g1 .+ t.*g2, q) for t in 0:q-1], [copy(g2)])
    for d in dirs
        p=findfirst(!=(0), d); p===nothing && continue
        dp=d[p]
        for val in -bnd:bnd
            val==0 && continue
            s=mod(val*inv_(dp), q)
            w=ctr(mod.(s.*d, q)); m=maximum(abs,w)
            if m<=bnd && m<bn && any(!=(0),w); bn=m; best=w; end
        end
    end
    best
end
function trial(n,seed; Mmul=1.2)
    rng=MersenneTwister(seed)
    f=Int[rand(rng,-5:5) for _ in 1:n]; f[1]=f[1]==0 ? 1 : f[1]
    fp=mod.(-sh(f,n÷2),q); k1=n÷2-1;k2=n-1;D=n*(n+1)÷2;M=Int(round(Mmul*D))
    meas=Vector{Int}[]
    for _ in 1:M
        call=rand(rng,1:2); k= call==1 ? k1 : k2; c=Int[rand(rng,0:q-1) for _ in 1:n]
        row=frow(c,k,n); s=mod(sum(row[j]*f[j] for j in 1:n),q); j0=findfirst(j->gcd(f[j],q)==1,1:n); p0=k-(j0-1)
        p0>=0 ? (c[p0+1]=mod(c[p0+1]-s*inv_(f[j0]),q)) : (c[p0+n+1]=mod(c[p0+n+1]+s*inv_(f[j0]),q))
        push!(meas, symmeas(frow(c,k1,n),n))    # solver only ever forms a_i (=row at k1)
    end
    A=reduce(vcat,[reshape(r,1,:) for r in meas])
    B=kernelbasis(A,D); length(B)==1 || return (false, length(B))
    S=tomat(B[1],n); cb=colbasis(S,n); length(cb)==2 || return (false,-2)
    w=shortvec(cb[1],cb[2],n)
    ok = w!==nothing && (w==f || w==Int.(.-f) || w==ctr(fp) || w==ctr(mod.(-fp,q)))
    (ok, 1)
end
function main()
    println("# blind_solve.jl -- recover f without which-call labels, via a rank-2 lift")
    @printf("%-4s %-6s %-8s %s\n","n","D","M","result (8 trials)")
    for n in (8,16,32)
        D=n*(n+1)÷2; M=Int(round(1.2*D)); ok=0
        for s in 1:8; r,_=trial(n,s; Mmul=1.4); r && (ok+=1); end
        @printf("%-4d %-6d %-8d recovered f blindly %d/8  (M = 1.4*n(n+1)/2 events)\n", n, D, M, ok)
    end
    println("\n# labelled channel needs O(n) events (section 6.5); blind needs O(n^2).")
    println("# at n=512 that is ~1.3e7 vs ~3.3e9 first-two events -- an n-fold overhead,")
    println("# not an impossibility.  The which-call disjunction is solved, not assumed away.")
end
main()
