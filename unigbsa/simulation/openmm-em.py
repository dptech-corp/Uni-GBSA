import os
import argparse

try:
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    from openmm.app.internal.pdbstructure import PdbStructure
    from openmm.app.forcefield import NonbondedGenerator
    from openmm import version
except:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    from simtk.openmm import unit
    from simtk.openmm.app.internal.pdbstructure import PdbStructure
    from simtk.openmm.app.forcefield import NonbondedGenerator
    from simtk.openmm import version
# Support Cythonized functions in OpenMM 7.3
# and also implementations in older versions.
if float(version.short_version[:3])<7.3:
    matchResidue = app.forcefield._matchResidue
else:
    try:
        from openmm.app.internal import compiled
    except ImportError:
        from simtk.openmm.app.internal import compiled
    matchResidue = compiled.matchResidueToTemplate

def energy_minimization(pdbfile, outfile=None, seed=None):
    file = open(pdbfile, 'r')
    structure = PdbStructure(file)
    pdb = app.PDBFile(structure)
    topology = pdb.topology
    positions = pdb.positions
    file.close()

    if outfile is None:
        outfile = os.path.split(pdbfile)[-1][:-4]+'_em.pdb'
    forcefield = createForceField(topology, False)
    system = forcefield.createSystem(topology)

    # If any heavy atoms were omitted, add them back to avoid steric clashes.

    #nonbonded = [f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce)]#[0]
    atomPositions = []
    for atom in topology.atoms():
        atomPositions.append(positions[atom.index])

    # For efficiency, only compute interactions that involve a new atom.

    #nonbonded.addInteractionGroup([atom.index for atom in topology.atoms()], range(system.getNumParticles()))

    # Do an energy minimization.

    integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 50*unit.femtosecond)
    if seed is not None:
        integrator.setRandomNumberSeed(seed)
    context = mm.Context(system, integrator)
    context.setPositions(atomPositions)
    mm.LocalEnergyMinimizer.minimize(context, maxIterations=10000)
    #s = app.simulation.Simulation(context)
    #s.step(1000)
    state = context.getState(getPositions=True)
    positions= state.getPositions()
    with open(outfile, 'w') as f:
        f.write("REMARK   1 energy minimization FROM: %s\n" % pdbfile)
        app.PDBFile.writeFile(topology, positions, f, True)

def createForceField(newTopology, water=False):
    """Create a force field to use for optimizing the positions of newly added atoms."""
    if water:
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        nonbonded = [f for f in forcefield._forces if isinstance(f, NonbondedGenerator)][0]
        radii = {'H':0.198, 'Li':0.203, 'C':0.340, 'N':0.325, 'O':0.299, 'F':0.312, 'Na':0.333, 'Mg':0.141,
                 'P':0.374, 'S':0.356, 'Cl':0.347, 'K':0.474, 'Br':0.396, 'Rb':0.527, 'I':0.419, 'Cs':0.605}
    else:
        forcefield = app.ForceField(os.path.join(os.path.dirname(__file__), '../data/soft.xml'))
    # The Topology may contain residues for which the ForceField does not have a template.
    # If so, we need to create new templates for them.
    atomTypes = {}
    bondedToAtom = []
    for atom in newTopology.atoms():
        bondedToAtom.append(set())
    for atom1, atom2 in newTopology.bonds():
        bondedToAtom[atom1.index].add(atom2.index)
        bondedToAtom[atom2.index].add(atom1.index)
    for residue in newTopology.residues():
        # Make sure the ForceField has a template for this residue.
        signature = app.forcefield._createResidueSignature([atom.element for atom in residue.atoms()])
        if signature in forcefield._templateSignatures:
            if any(matchResidue(residue, t, bondedToAtom) is not None for t in forcefield._templateSignatures[signature]):
                continue
        # Create a new template.
        resName = "extra_"+residue.name
        template = app.ForceField._TemplateData(resName)
        forcefield._templates[resName] = template
        indexInResidue = {}
        for atom in residue.atoms():
            element = atom.element
            typeName = 'extra_'+element.symbol
            if element not in atomTypes:
                atomTypes[element] = app.ForceField._AtomType(typeName, '', 0.0, element)
                forcefield._atomTypes[typeName] = atomTypes[element]
                if water:
                    # Select a reasonable vdW radius for this atom type.
                    if element.symbol in radii:
                        sigma = radii[element.symbol]
                    else:
                        sigma = 0.5
                    nonbonded.registerAtom({'type':typeName, 'charge':'0', 'sigma':str(sigma), 'epsilon':'0'})
            indexInResidue[atom.index] = len(template.atoms)
            template.atoms.append(app.ForceField._TemplateAtomData(atom.name, typeName, element))
        for atom in residue.atoms():
            for bondedTo in bondedToAtom[atom.index]:
                if bondedTo in indexInResidue:
                    b = (indexInResidue[atom.index], indexInResidue[bondedTo])
                    if b[0] < b[1]:
                        template.bonds.append(b)
                        template.atoms[b[0]].bondedTo.append(b[1])
                        template.atoms[b[1]].bondedTo.append(b[0])
                else:
                    b = indexInResidue[atom.index]
                    template.externalBonds.append(b)
                    template.atoms[b].externalBonds += 1
        if signature in forcefield._templateSignatures:
            forcefield._templateSignatures[signature].append(template)
        else:
            forcefield._templateSignatures[signature] = [template]
    return forcefield

def will_restrain(atom, rset):
  """Returns True if the atom will be restrained by the given restraint set."""

  if rset == "heavy":
    return atom.element.name != "hydrogen"
  elif rset == "ca":
    return atom.name == "CA"

def minimize(pdbfile, outfile, gmxtop=False, use_gpu=False, rset='ca', exclude_residues={}):

    pdb = app.PDBFile(pdbfile)
    if gmxtop:
        ff = app.GromacsTopFile(gmxtop)
        system = ff.createSystem(nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)
    else:
        ff = app.ForceField('amber99sb.xml', 'tip3p.xml')
        system = ff.createSystem(pdb.topology,
                             nonbondedMethod=app.CutoffNonPeriodic,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds)
    integrator = mm.LangevinMiddleIntegrator(300.0 * unit.kelvin,
                                        5.0 / unit.picosecond,
                                        1.0 * unit.femtosecond)
    platform = mm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
    # integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 2*unit.femtosecond)

    # harmonic restrain
    if rset:
        """Adds a harmonic potential that restrains the system to a structure."""
        assert rset in ["heavy", "ca"]
        restraint = mm.CustomExternalForce('0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2')
        system.addForce(restraint)
        restraint.addGlobalParameter('k', 10.0*unit.kilocalorie_per_mole/unit.angstrom**2)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')

        for atom in pdb.topology.atoms():
            if atom.residue.chain.id in exclude_residues:
                key = '%s%s'%(atom.residue.id, atom.residue.insertionCode.strip())
                if key in exclude_residues[atom.residue.chain.id]:
                    print('Exclude atom:', atom)
                    continue
            if will_restrain(atom, rset):
                restraint.addParticle(atom.index, pdb.positions[atom.index])

    context = mm.Context(system, integrator, platform)
    context.setPositions(pdb.getPositions())
    mm.LocalEnergyMinimizer_minimize(context, tolerance=2.39*unit.kilocalorie_per_mole, maxIterations=0)
    state = context.getState(getEnergy=True, getPositions=True, groups={0})
    with open(outfile, "w") as f:
        app.PDBFile.writeFile(pdb.topology,
                              state.getPositions(asNumpy=True),
                              f,
                              keepIds=True)

def main():
    parser = argparse.ArgumentParser(description='Do energy minimization for input pdbfile.')
    parser.add_argument('-i', dest='inp', help='Input PDB file.', required=True)
    parser.add_argument('-o', dest='out', help='Output PDB file.', default=None)
    parser.add_argument('-p', dest='top', help='Gromacs topology file.', default=None)
    parser.add_argument('-restrain', help='Atoms to restrain. default: none', choices=['ca', 'heavy', 'none'], default=None)
    parser.add_argument('-exclude_residues', help='Residues exclude for restrain. format: chainResid[Insercode], eg: A12,B13A', nargs='+', default=[])
    args = parser.parse_args()
    pdbfile, outfile, rset, top = args.inp, args.out, args.restrain, args.top
    excludeResidues = {}
    for item in args.exclude_residues:
        chain = item[0]
        if chain in excludeResidues:
            excludeResidues[chain].append(item[1:])
        else:
            excludeResidues[chain] = [item[1:]]

    #energy_minimization(pdbfile, outfile)
    if rset:
        rset = rset.lower()
    minimize(pdbfile, outfile, gmxtop=top, use_gpu=False, rset=rset, exclude_residues=excludeResidues)

if __name__ == '__main__':
    main()
