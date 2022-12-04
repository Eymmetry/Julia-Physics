mutable struct Ball
    v :: Float64
    mass :: Float64
end

function check_collision!(small::Ball, large::Ball)
    r = large.mass / small.mass
    v = small.v
    V = large.v
    if v > V
        V_next = (2v+(r-1)V)/(r+1)
        v_next = ((1-r)v+2r*V)/(r+1)
        collision = true

        small.v = v_next
        large.v = V_next
        return collision
    elseif  v < 0
        collision = true
        small.v = -v
        return collision
    else
        collision =false
        return collision
    end文献
end

function ballcollison(N)
    v0 = 0.0
    V0 = -1.0
    m = 1
    M = m*100^N
    smallball = Ball(v0,m)
    largeball = Ball(V0,M)
    collision = true
    cnt = 0
    while collision
        collision = check_collision!(smallball,largeball)
        if collision
            cnt += 1
        end
    end
    return cnt/10^N
end

N = 6
p = ballcollison(N)
println(p)
