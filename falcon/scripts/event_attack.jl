#!/usr/bin/env julia
#
# event_attack.jl -- the first-two event channel, END TO END, on REAL hashed
# syndromes: recover the key from integer-centre events and forge.
#
#     julia --project=falcon falcon/scripts/event_attack.jl [n] [scan_cap]
#
# scripts/event_solve.jl established the solve but drew syndromes UNIFORMLY --
# it did not hash messages.  That was its stated idealization (section 6.5).
# This script removes it: every syndrome is a real hash_to_point(message, salt)
# output, exactly as a signer produces, and the events are the real condition
#
#     (c*f)[n/2-1] = 0   (call 1)   or   (c*f)[n-1] = 0   (call 2)   (mod q),
#
# which section 6.4 proved is the exact first-two sampler-centre integer
# condition.  Both are linear equations on f with public coefficients, so a
# mix of call-1 and call-2 events STACKS into one F_q system.  We collect n-1
# events, solve, recover f, complete the trapdoor basis (solving fG - gF = q,
# the same descent key generation runs), sign a NEW message, and verify it
# against the victim's public key.
#
# What is idealised here and what is not.  The syndromes are real hashes (the
# section 6.5 idealisation is removed).  The event LABEL -- which of the two
# calls straddled -- is taken as known here; scripts/window.jl and section 6.7
# measure the cross-build ORACLE that detects a first-two divergence, and the
# residual for a fully blind attacker is disambiguating call 1 from call 2
# (a one-of-two disjunction; naive iterative solving does not crack it, so we
# report it as an open problem rather than claim it).
#
# CONSTANT TIME: irrelevant; attacker side.

using Falcon, Printf, Random
const F = Falcon

function form_row(c::Vector{Int}, k::Int, n::Int, q::Int)
    row = Vector{Int}(undef, n)
    @inbounds for j in 0:(n-1)
        idx = k - j
        row[j+1] = idx >= 0 ? mod(c[idx+1], q) : mod(-c[idx+n+1], q)
    end
    return row
end
function rref_modq(rows, n, q)
    A=[copy(r) for r in rows]; piv=Int[]; r=1
    for col in 1:n
        p=findfirst(i->A[i][col]%q!=0, r:length(A)); p===nothing&&continue; p+=r-1
        A[r],A[p]=A[p],A[r]; iv=invmod(A[r][col],q); A[r]=Int[mod(x*iv,q) for x in A[r]]
        for i in 1:length(A); i==r&&continue; fc=A[i][col]; fc==0&&continue
            A[i]=Int[mod(A[i][t]-fc*A[r][t],q) for t in 1:n]; end
        push!(piv,col); r+=1; r>length(A)&&break
    end
    r-1,A,piv
end
function kernel_vec(A,piv,n,q)
    free=setdiff(1:n,piv); length(free)==1||return nothing; fc=free[1]
    v=zeros(Int,n); v[fc]=1
    for (r,c) in enumerate(piv); v[c]=mod(-A[r][fc],q); end; v
end
ctr(x,q)=(y=mod(x,q); y>q÷2 ? y-q : y)

function main()
    n = length(ARGS)>=1 ? parse(Int,ARGS[1]) : 512
    cap = length(ARGS)>=2 ? parse(Int,ARGS[2]) : 40_000_000
    p = n==512 ? FALCON_512 : FALCON_1024
    q=p.q; k1=n÷2-1; k2=n-1
    rb = chacha20(collect(UInt8, 0x00:0x37))
    sk = falcon_keygen(n, kk->randombytes!(rb,kk))[1]
    f = Int.(sk.f); g = Int.(sk.g)
    h = F.polydivq(Int[mod(c,q) for c in g], Int[mod(c,q) for c in f])
    @printf("# event_attack.jl  n=%d  q=%d   target max|f|=%d\n", n, q, maximum(abs,f))
    @printf("# scanning real hash_to_point syndromes for first-two events\n\n")

    rows=Vector{Int}[]; tried=0; c1=0; c2=0
    rank=0; A=nothing; piv=nothing
    while rank < n-1 && tried < cap
        tried += 1
        msg = collect(codeunits("m/$tried"))
        salt = shake256(collect(codeunits("s/$tried")), SALT_LEN)
        c = hash_to_point(msg, salt, n; q=q)
        r1 = form_row(c,k1,n,q); v1 = mod(sum(r1[j]*f[j] for j in 1:n), q)
        r2 = form_row(c,k2,n,q); v2 = mod(sum(r2[j]*f[j] for j in 1:n), q)
        added=false
        if v1==0; push!(rows,r1); c1+=1; added=true; end
        if v2==0; push!(rows,r2); c2+=1; added=true; end
        if added && (length(rows)%128==0 || length(rows)>=n-1)
            rank,A,piv = rref_modq(rows,n,q)
        end
    end
    @printf("scanned %d messages; %d call-1 events + %d call-2 events = %d rows\n",
            tried, c1, c2, length(rows))
    @printf("observed event rate %.3g  (predicted 2/q = %.3g)\n", (c1+c2)/tried, 2/q)
    rank,A,piv = rref_modq(rows,n,q)
    @printf("system rank %d of %d\n", rank, n)
    v = kernel_vec(A,piv,n,q)
    v===nothing && (println("kernel not 1-dim; scan more"); return)
    best=nothing; bm=typemax(Int)
    for s in 1:q-1
        w=Int[ctr(s*v[i],q) for i in 1:n]; m=maximum(abs,w)
        m<bm && (bm=m; best=w); bm<=maximum(abs,f)&&break
    end
    ok = best==f || best==.-f
    @printf("\nrecovered vector max|.| = %d ; matches secret f: %s%s\n",
            bm, ok, best==.-f ? " (up to sign)" : "")
    ok || (println("recovery mismatch"); return)

    # forge: complete the basis, sign a fresh message, verify against the pubkey
    fr = best==f ? f : Int.(.-f); gg = best==f ? g : Int.(.-g)
    try
        bF,bG = F.ntru_solve(BigInt.(fr), BigInt.(gg); q=q)
        fsk = expand_privkey(fr, gg, Int.(bF), Int.(bG), p)
        pk = FalconPublicKey(p, h)
        m = collect(codeunits("forgery via first-two events"))
        rr = chacha20(shake256(collect(codeunits("forge-fte")),56))
        sig = falcon_sign(fsk, m, x->randombytes!(rr,x))
        acc = falcon_verify(pk, m, sig)
        @printf("forged a NEW message with the recovered key: %s\n",
                acc ? "*** the victim's public key ACCEPTS it ***" : "rejected")
    catch e
        @printf("basis completion failed: %s\n", sprint(showerror,e))
    end
end
main()
