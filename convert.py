import glob
import os, sys
from CHARMM import CHARMM

inputpsf = sys.argv[1]
inputpdb = sys.argv[2]

log = open('out', 'w')
charmm = CHARMM(output=log, exe='/Users/sunhwan/local/charmm/c40a1')
charmm.sendCommand("stream toppar.str")
charmm.sendCommand("stream toppar/toppar_water_ions.str")

charmm.sendCommand("""
    bomlev -5
    ioformat noext
    read psf card name %(inputpsf)s xplor
    read coor pdb name %(inputpdb)s resid
    ioformat exte
    """ % locals())

# join water segments
charmm.sendCommand("""
    define jnk sele resn tip3 end
    set nsel = ?nsel
    set firstsegid = ?selsegi
    coor stat sele segid @firstsegid end
    set firstnsel = ?nsel""")

natom = charmm.params['nsel']
nsel = charmm.params['firstnsel'] + 1
while 1:
    charmm.sendCommand("""
        define jnk sele resn tip3 .subset. %(nsel)d end
        set segid = ?selsegi
        coor stat sele segid @segid end
        calc nsel = %(nsel)d + ?nsel
        join @firstsegid @segid renumber""" % locals())
    nsel = charmm.params['nsel']
    if nsel > natom: break

# rename tip3 to tip4
charmm.sendCommand("""
    rename resname tip4 sele segid @firstsegid end
    """)
water_segid = charmm.params['firstsegid']

# split individual segments
charmm.sendCommand("""
    define jnk sele all end
    set nsel = ?nsel""")
natom = charmm.params['nsel']

segids = []
nsel = 1
while 1:
    charmm.sendCommand("""
        define jnk sele all .subset. %(nsel)d end
        set segid = ?selsegi
        coor stat sele segid @segid end
        calc nsel = %(nsel)d + ?nsel
        write coor card name /tmp/@segid.crd sele segid @segid end
        * @segid
        *
        """ % locals())

    segids.append(charmm.params['segid'])
    nsel = charmm.params['nsel']
    if nsel > natom: break


# read each segments back
charmm = CHARMM(output=log, exe='/Users/sunhwan/local/charmm/c40a1')
charmm.sendCommand("stream toppar.str")
charmm.sendCommand("stream toppar/toppar_water_ions_tip4p.str")
 
psffile = "%s_tip4.psf" % '.'.join(os.path.basename(inputpsf).split('.')[:-1])
pdbfile = "%s_tip4.pdb" % '.'.join(os.path.basename(inputpdb).split('.')[:-1])

for segid in segids:
    if segid == water_segid:
        first = 'None'
        last = 'None'
        extra = 'noangle nodihedral'
    elif segid == 'SOD':
        first = 'None'
        last = 'None'
        extra = 'noangle nodihedral'
    else:
        first = 'None'
        last = 'None'
        extra = ''

    charmm.sendCommand("""
        read sequ coor name /tmp/%(segid)s.crd
        generate %(segid)s first %(first)s last %(last)s setup warn %(extra)s

        read coor card name /tmp/%(segid)s.crd
        """ % locals())

    if segid == water_segid:
        # build OM position
        #charmm.sendCommand("stream build_tip4_om.str %(segid)s" % locals())
        charmm.sendCommand("""
            set segid = %(segid)s
            define tip4 sele type oh2 .and. segid %(segid)s end
            set nwater = ?nsel
            """ % locals())
        nwater = int(charmm.params["nwater"])

        for i in range(nwater):
            charmm.sendCommand("""
               calc icnt = %(i)d + 1
               coor stat sele type oh2 .and. resid @icnt .and. segid @segid end
               set xoh2 = ?xave
               set yoh2 = ?yave
               set zoh2 = ?zave

               coor stat sele type h* .and. resid @icnt .and. segid @segid end
               set xh2 = ?xave
               set yh2 = ?yave
               set zh2 = ?zave   

               calc bisectorx = @xoh2 - @xh2
               calc bisectory = @yoh2 - @yh2
               calc bisectorz = @zoh2 - @zh2
               calc r = sqrt(@bisectorx * @bisectorx + @bisectory * @bisectory + @bisectorz * @bisectorz)

               calc omx = @xoh2 + @bisectorx / @r * 0.15
               calc omy = @yoh2 + @bisectory / @r * 0.15
               calc omz = @zoh2 + @bisectorz / @r * 0.15

               coor set xdir @omx ydir @omy zdir @omz sele type om .and. resid @icnt .and. segid @segid end
            """ % locals())

        # rename water segid
        charmm.sendCommand("""
            rename segid TIP4 sele segid %(segid)s end
            """ % locals())

charmm.sendCommand("""
    write psf card name %(psffile)s xplor
    * tip4
    *
    """ % locals())

charmm.sendCommand("""
    write coor pdb name %(pdbfile)s
    * tip4
    *
    """ % locals())
