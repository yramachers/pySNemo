import pysnemo.utility.euclid as EU
from pysnemo.utility.interval import Interval

# consider lines only for alternative geometries
def get_calorimeter_plane(info):
    type = info[1] # calo type
    side = info[3] # the tracker side
    wall = info[6] # the wall
    
    # Main Wall
    if type==0:
        # main calo sitting at x = +- 27.5738cm
        if side == 0:
            calo = EU.Plane(EU.Point3(-275.738,0.0,0.0),EU.Point3(-275.738,1.0,0.0),EU.Point3(-275.738,0.0,1.0)) 
        else:
            calo = EU.Plane(EU.Point3(275.738,0.0,0.0),EU.Point3(275.738,1.0,0.0),EU.Point3(275.738,0.0,1.0)) 

    # X Wall
    elif type==1:
        # xwall calo sitting at y = +- 250.55cm
        if wall == 0:
            calo = EU.Plane(EU.Point3(0.0, -2505.5,0.0),EU.Point3(1.0, -2505.5,0.0),EU.Point3(0.0,-2505.5,1.0)) 
        else:
            calo = EU.Plane(EU.Point3(0.0, 2505.5,0.0),EU.Point3(1.0, 2505.5,0.0),EU.Point3(0.0,2505.5,1.0)) 

    # Gamma Veto
    elif type==2:
        # gveto calo sitting at z = +- 155.0cm
        if wall == 0:
            calo = EU.Plane(EU.Point3(1.0,0.0, -1550.0),EU.Point3(0.0,1.0, -1550.0),EU.Point3(0.0,0.0, -1550.0)) 
        else:
            calo = EU.Plane(EU.Point3(1.0,0.0, 1550.0),EU.Point3(0.0,1.0, 1550.0),EU.Point3(0.0,0.0, 1550.0)) 

    return calo


def get_ellipse(interpoints, type):
    ell = []
    # Main Wall
    if type==0:
        ell.append(interpoints[0].y)
        ell.append(interpoints[0].z)
        ell.append(interpoints[1].y)
        ell.append(interpoints[2].y)
        ell.append(interpoints[3].z)
        ell.append(interpoints[4].z)

    # X Wall
    elif type==1:
        ell.append(interpoints[0].x)
        ell.append(interpoints[0].z)
        ell.append(interpoints[1].x)
        ell.append(interpoints[2].x)
        ell.append(interpoints[3].z)
        ell.append(interpoints[4].z)

    # Gamma Veto
    elif type==2:
        ell.append(interpoints[0].x)
        ell.append(interpoints[0].y)
        ell.append(interpoints[1].x)
        ell.append(interpoints[2].x)
        ell.append(interpoints[3].y)
        ell.append(interpoints[4].y)

    return ell



def calopar(lines, info):
    caloplane = get_calorimeter_plane(info)

    interpoints = [] # center, left, right, top, bottom
    for line in lines:
        interpoints.append(line.intersect(caloplane)) # 5 Point3 objects

    ellipse = get_ellipse(interpoints, info[1])
    return ellipse


def foilpar(linepar):
    best = linepar[0]
    err = linepar[1]

    icxy = best[0]
    icxz = best[2]
    err_icxy = err[0]
    err_icxz = err[2]

    ellipse = [icxy,icxz,icxy+err_icxy,icxy-err_icxy,icxz+err_icxz,icxz-err_icxz]
    return ellipse


def create_lines(linepar):
    '''
    Expect 4 line parameters, i.e. icxy, slxy, icxz, slxz
    to create a line3 object used in extrapolation.
    Intercepts are taken at foil x=0 by construction.
    '''
    best = linepar[0]
    err = linepar[1]

    icxy = best[0]
    slxy = best[1]
    icxz = best[2]
    slxz = best[3]
    err_icxy = err[0]
    err_slxy = err[1]
    err_icxz = err[2]
    err_slxz = err[3]

    # center
    p = EU.Point3(0.0,icxy,icxz)
    vec = EU.Vector3(1.0,slxy,slxz)
    lc = EU.Line3(p,vec)
    # left
    if (slxy>=0.0):
        p = EU.Point3(0.0,icxy-err_icxy,icxz)
        vec = EU.Vector3(1.0,slxy+err_slxy,slxz)
        ll = EU.Line3(p,vec)
    else:
        p = EU.Point3(0.0,icxy-err_icxy,icxz)
        vec = EU.Vector3(1.0,slxy-err_slxy,slxz)
        ll = EU.Line3(p,vec)        
    # right
    p = EU.Point3(0.0,icxy+err_icxy,icxz)
    vec = EU.Vector3(1.0,slxy-err_slxy,slxz)
    lr = EU.Line3(p,vec)
    # top
    if (slxz>=0.0):
        p = EU.Point3(0.0,icxy,icxz-err_icxz)
        vec = EU.Vector3(1.0,slxy,slxz+err_slxz)
        lt = EU.Line3(p,vec)
    else:
        p = EU.Point3(0.0,icxy,icxz-err_icxz)
        vec = EU.Vector3(1.0,slxy,slxz-err_slxz)
        lt = EU.Line3(p,vec)
    # bottom
    p = EU.Point3(0.0,icxy,icxz+err_icxz)
    vec = EU.Vector3(1.0,slxy,slxz-err_slxz)
    lb = EU.Line3(p,vec)
    return [lc,ll,lr,lt,lb]



def hit_the_bar(intA, intB, info):
#    print 'Got: type=%d, side=%d, wall=%d, column=%d, row=%d'%(info[0],info[1],info[4],info[2],info[3])
    type = info[0]
    side = info[1]
    column = info[2]
    row = info[3]
    wall = info[4]
    
    if type==0:
        yinit = column * 67.0 - 2546.0 # 65 for modules + 1mm offset
        #zinit = row * 259.0 - 1683.5 # 259 for modules + offset
        
        dy = Interval(yinit + 0.99, yinit + 1.01 + 65.0)
        dz = Interval(-1683.5, 1683.5) # one long bar calo
#        print 'Got intervals y, z: '
#        print intA
#        print dy
#        print intB
#        print dz
        
        
        if (intA.overlap(dy)) and (intB.overlap(dz)):
            return True
        else:
            return False
        

    elif type==1:
        if side == 0:
            xinit = -(column * 121.0 + 29.) # 120 for module + 0.5mm gap + offset
            dx = Interval(xinit - 120.0 - 0.51, xinit - 0.49)
        else:
            xinit = column * 121.0 + 29. # 120 for module + 0.5mm gap + offset
            dx = Interval(xinit + 0.49, xinit + 0.51 + 120.0)
            
        zinit = row * 212.0 - 1696.0 # 212 for modules + offset

        dz = Interval(zinit + 1.749, zinit + 1.75 + 208.51)
#        print 'Got intervals x, z: '
#        print intA
#        print dx
#        print intB
#        print dz


        if (intA.overlap(dx)) and (intB.overlap(dz)):
            return True
        else:
            return False
		
    else:
        leftgv = Interval(0,8) # 8 not included!
        if side == 0:
            xinit = - 11.9 # 11.9mm offset to source
            dx = Interval(xinit - 1.01 - 240.0, xinit - 0.99)
        else:
            xinit = 11.9 # 11.9mm offset to source
            dx = Interval(xinit + 0.99, xinit + 1.01 + 240.0)
            
        if column in leftgv: # in first block of 8
            yinit = column * 311.5 - 2497.25 # 311.5 for modules + offset
        else:
            yinit = column * 311.5 - 2497.25 + 10.5# 311.5 for modules + offset
            
        dy = Interval(yinit + 1.749, yinit + 1.75 + 308.01)
#        print 'Got intervals x, y: '
#        print dx
#        print dy
		
        if (intA.overlap(dx)) and (intB.overlap(dy)):
            return True
        else:
            return False
        


def particle_test(calopar, calo_mi):
    # geometry first
    # the calo hit
    calo_type = calo_mi[1]
    calo_side = calo_mi[3]
    calo_column = calo_mi[4]
    calo_row = calo_mi[5]
    calo_wall = calo_mi[6]

    info = [calo_type, calo_side, calo_column, calo_row, calo_wall]

    # the suspect particle hit
    if calopar[2]<=calopar[3]:
        particlehitA = Interval(calopar[2],calopar[3])
    else:
        particlehitA = Interval(calopar[3],calopar[2])
    if calopar[4]<=calopar[5]:
        particlehitB = Interval(calopar[4],calopar[5])
    else:
        particlehitB = Interval(calopar[5],calopar[4])
        
    if hit_the_bar(particlehitA, particlehitB, info):
        return True
    else:
        return False
    
    


def extrapolate_line(lines,fitside,calohit):
    '''
    Input one line in list to extrapolate to left and right
    '''
    calo_mi = calohit.meta_info
    if not fitside==calo_mi[3]: # calo and fit on opposite tracker halfs
        return []
    
    linepar = lines[0]
    ll = create_lines(linepar)
    f = foilpar(linepar)
    c = calopar(ll, calo_mi)
    
    if particle_test(c,calo_mi): # is a particle
        #            print 'Is a particle.'
        #            print f, c
        return [f, c] # list of tuples for foil and calo
    else:
        #            print 'Not a particle.'
        return [] # was path leading nowhere, no associated calo hit
    
