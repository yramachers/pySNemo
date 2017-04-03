import pysnemo.utility.euclid as EU
import pysnemo.utility.helix as HL
from pysnemo.utility.interval import Interval

def get_calorimeter_plane(info):
    type = info[1] # calo type
    side = info[3] # the tracker side
    wall = info[6] # the wall
    
    # Main Wall
    if type==0:
        # main calo sitting at x = +- 43.5cm
        if side == 0:
            calo = EU.Plane(EU.Point3(-435,0.0,0.0),EU.Point3(-435,1.0,0.0),EU.Point3(-435,0.0,1.0)) 
        else:
            calo = EU.Plane(EU.Point3(435,0.0,0.0),EU.Point3(435,1.0,0.0),EU.Point3(435,0.0,1.0)) 

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


def helix_calopar(hfit, info):
    caloplane = get_calorimeter_plane(info)
    type = info[1]

    helices = []
    helix = hfit.hel # best fit helix
    helices.append(helix)
    par = hfit.par
    err = hfit.errors
    # vary omega and tanlambda
    pl = (par[0],par[1],par[2],par[3]+err[3],par[4])
    helices.append(HL.Par5Helix(pl,hfit.bf))
    pr = (par[0],par[1],par[2],par[3]-err[3],par[4])
    helices.append(HL.Par5Helix(pr,hfit.bf))
    pt = (par[0],par[1],par[2],par[3],par[4]+err[4])
    helices.append(HL.Par5Helix(pt,hfit.bf))
    pb = (par[0],par[1],par[2],par[3],par[4]-err[4])
    helices.append(HL.Par5Helix(pb,hfit.bf))

    interpoints = []
    if type==0: # main
        ponp = caloplane._get_point()
        planetup = (ponp.x,ponp.y,1.0,0.0) # xaxis is normal
        for h in helices:
            tup = h.intersectionXY(planetup)
            if tup is None: # intersection failed
                return None
            p = EU.Point3(tup[0], tup[1], tup[2])
            interpoints.append(p)
    elif type==1: # xwall
        ponp = caloplane._get_point()
        planetup = (ponp.x,ponp.y,0.0,1.0) # yaxis is normal
        for h in helices:
            tup = h.intersectionXY(planetup)
            if tup is None: # intersection failed
                return None
            p = EU.Point3(tup[0], tup[1], tup[2])
            interpoints.append(p)
    else:
        zplane = caloplane._get_point().z
        for h in helices:
            tup = h.intersectionZ(zplane)
            if tup is None: # intersection failed
                return None
            p = EU.Point3(tup[0], tup[1], tup[2])
            interpoints.append(p)

    ellipse = get_ellipse(interpoints, type)
    return ellipse



def helix_foilpar(hfit):
    par = hfit.par
    err = hfit.errors
    # vary start point parameters, y, z +- err each
    ell = [par[1],par[2],par[1]+err[1],par[1]-err[1],par[2]+err[2],par[2]-err[2]]
    return ell



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



def hit_the_block(intA, intB, info):
#    print 'Got: type=%d, side=%d, wall=%d, column=%d, row=%d'%(info[0],info[1],info[4],info[2],info[3])
    type = info[0]
    side = info[1]
    column = info[2]
    row = info[3]
    wall = info[4]
    
    if type==0:
        yinit = column * 259.0 - 2590.0 # 259 for modules + offset
        zinit = row * 259.0 - 1683.5 # 259 for modules + offset
        
        dy = Interval(yinit + 1.49, yinit + 1.5 + 256.01)
        dz = Interval(zinit + 1.49, zinit + 1.5 + 256.01)
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
            xinit = -(column * 203.0 + 29.) # 202 for module + 1mm gap + offset
            dx = Interval(xinit - 200.0 - 1.01, xinit - 0.99)
        else:
            xinit = column * 203.0 + 29. # 202 for module + 1mm gap + offset
            dx = Interval(xinit + 0.99, xinit + 1.01 + 200.0)
            
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
            xinit = - 4.9 # 4.9mm offset to source
            dx = Interval(xinit - 1.01 - 290.0, xinit - 0.99)
        else:
            xinit = 4.9 # 4.9mm offset to source
            dx = Interval(xinit + 0.99, xinit + 1.01 + 290.0)
            
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
                                  
    if hit_the_block(particlehitA, particlehitB, info):
        return True
    else:
        return False
    
    


def extrapolate_helix(hfit,fitside,calohit):
    '''
    Input is HelixFit object for extrapolation

    '''
    calo_mi = calohit.meta_info
    if not fitside==calo_mi[3]: # calo and fit on opposite tracker halfs
        return []

    c = helix_calopar(hfit, calo_mi)
    if c is None: # intersection failed
        return []

    f = helix_foilpar(hfit)
#    print 'HELIX calo par: ', c

    if particle_test(c,calo_mi): # is a particle
        binflag = False # no kinks here
        angles = [] # empty here
        return [f, c, binflag, angles] # list of tuples for foil and calo
    else:
        return [] # was path leading nowhere, no associated calo hit



def extrapolate_line(lines,fitside,calohit):
    '''
    Input either one line in list to extrapolate to left and right
    or a broken line giving two line segments in list for extrapolation
    to left and right. 
    '''
    calo_mi = calohit.meta_info
    if not fitside==calo_mi[3]: # calo and fit on opposite tracker halfs
        return []
    
    if len(lines)<2:
        linepar = lines[0]
        ll = create_lines(linepar)
        f = foilpar(linepar)
        c = calopar(ll, calo_mi)

        if particle_test(c,calo_mi): # is a particle
            binflag = False # no kinks here
            angles = [] # empty here
#            print 'Is a particle.'
#            print f, c
            return [f, c, binflag, angles] # list of tuples for foil and calo
        else:
#            print 'Not a particle.'
            return [] # was path leading nowhere, no associated calo hit

    else: # this is the broken line fit result, no Migrad
        linepar = lines[0]
        f = foilpar(linepar) # extrapolator to foil
        linepar = lines[1]
        ll = create_lines(linepar)
        c = calopar(ll, calo_mi) # extrapolator to calo

        if particle_test(c,calo_mi): # is a particle
            binflag = True
            angles = linepar[-1] # any linepar will do, angles are stored in all
            return [f, c, binflag, angles]
        else:
            return []
