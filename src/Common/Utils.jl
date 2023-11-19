function cart2sph(x, y, z)
    azimuth = atan(y, x)
    elevation = atan(z, sqrt(x^2 + y^2))
    r = sqrt(x^2 + y^2 + z^2)

    return (azimuth, elevation, r)
end

function cart2sph(pos::AbstractArray)
    azimuth = atan(pos[2], pos[1])
    elevation = atan(pos[3], sqrt(pos[1]^2 + pos[2]^2))
    r = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)

    return [azimuth, elevation, r]
end

function sph2cart(az, el, r)
    x = r * cos(el) * cos(az)
    y = r * cos(el) * sin(az)
    z = r * sin(el)

    return (x, y, z)
end

function sph2cart(sph::AbstractArray)
    x = sph[3] * cos(sph[2]) * cos(sph[1])
    y = sph[3] * cos(sph[2]) * sin(sph[1])
    z = sph[3] * sin(sph[2])

    return [x, y, z]
end

function angVec(a::AbstractVector, b::AbstractVector)
    return acos(clamp(dot(a, b) / (norm(a) * norm(b)), -1, 1))
end