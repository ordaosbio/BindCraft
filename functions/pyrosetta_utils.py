####################################
################ PyRosetta functions
####################################
### Import dependencies
import os
#import pyrosetta as pr
#from pyrosetta.rosetta.core.kinematics import MoveMap
#from pyrosetta.rosetta.core.select.residue_selector import ChainSelector
#from pyrosetta.rosetta.protocols.simple_moves import AlignChainMover
#from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
#from pyrosetta.rosetta.protocols.relax import FastRelax
#from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
#from pyrosetta.rosetta.core.select import get_residues_from_subset
#from pyrosetta.rosetta.core.io import pose_from_pose
#from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from .generic_utils import clean_pdb
from .biopython_utils import hotspot_residues
from biotite.structure import superimpose_homologs, rmspd
from biotite.structure.io import load_structure, save_structure
from alphafold.relax import relax
from alphafold.common import protein, residue_constants

# Rosetta interface scores
def score_interface(pdb_file, binder_chain="B"):
    # load pose
    #pose = pr.pose_from_pdb(pdb_file)

    # analyze interface statistics
    #iam = InterfaceAnalyzerMover()
    #iam.set_interface("A_B")
    #scorefxn = pr.get_fa_scorefxn()
    #iam.set_scorefunction(scorefxn)
    #iam.set_compute_packstat(True)
    #iam.set_compute_interface_energy(True)
    #iam.set_calc_dSASA(True)
    #iam.set_calc_hbond_sasaE(True)
    #iam.set_compute_interface_sc(True)
    #iam.set_pack_separated(True)
    #iam.apply(pose)

    # Initialize dictionary with all amino acids
    interface_AA = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}

    # Initialize list to store PDB residue IDs at the interface
    interface_residues_set = hotspot_residues(pdb_file, binder_chain)
    interface_residues_pdb_ids = []

    # Iterate over the interface residues
    for pdb_res_num, aa_type in interface_residues_set.items():
        # Increase the count for this amino acid type
        interface_AA[aa_type] += 1

        # Append the binder_chain and the PDB residue number to the list
        interface_residues_pdb_ids.append(f"{binder_chain}{pdb_res_num}")

    # count interface residues
    interface_nres = len(interface_residues_pdb_ids)

    # Convert the list into a comma-separated string
    interface_residues_pdb_ids_str = ','.join(interface_residues_pdb_ids)

    # Calculate the percentage of hydrophobic residues at the interface of the binder
    hydrophobic_aa = set('ACFILMPVWY')
    hydrophobic_count = sum(interface_AA[aa] for aa in hydrophobic_aa)
    if interface_nres != 0:
        interface_hydrophobicity = (hydrophobic_count / interface_nres) * 100
    else:
        interface_hydrophobicity = 0

    # retrieve statistics
    #interfacescore = iam.get_all_data()
    #interface_sc = interfacescore.sc_value # shape complementarity
    #interface_interface_hbonds = interfacescore.interface_hbonds # number of interface H-bonds
    #interface_dG = iam.get_interface_dG() # interface dG
    #interface_dSASA = iam.get_interface_delta_sasa() # interface dSASA (interface surface area)
    #interface_packstat = iam.get_interface_packstat() # interface pack stat score
    #interface_dG_SASA_ratio = interfacescore.dG_dSASA_ratio * 100 # ratio of dG/dSASA (normalised energy for interface area size)
    #buns_filter = XmlObjects.static_get_filter('<BuriedUnsatHbonds report_all_heavy_atom_unsats="true" scorefxn="scorefxn" ignore_surface_res="false" use_ddG_style="true" dalphaball_sasa="1" probe_radius="1.1" burial_cutoff_apo="0.2" confidence="0" />')
    #interface_delta_unsat_hbonds = buns_filter.report_sm(pose)

    #if interface_nres != 0:
    #    interface_hbond_percentage = (interface_interface_hbonds / interface_nres) * 100 # Hbonds per interface size percentage
    #    interface_bunsch_percentage = (interface_delta_unsat_hbonds / interface_nres) * 100 # Unsaturated H-bonds per percentage
    #else:
    #    interface_hbond_percentage = None
    #    interface_bunsch_percentage = None

    # calculate binder energy score
    #chain_design = ChainSelector(binder_chain)
    #tem = pr.rosetta.core.simple_metrics.metrics.TotalEnergyMetric()
    #tem.set_scorefunction(scorefxn)
    #tem.set_residue_selector(chain_design)
    #binder_score = tem.calculate(pose)

    # calculate binder SASA fraction
    #bsasa = pr.rosetta.core.simple_metrics.metrics.SasaMetric()
    #bsasa.set_residue_selector(chain_design)
    #binder_sasa = bsasa.calculate(pose)

    #if binder_sasa > 0:
    #    interface_binder_fraction = (interface_dSASA / binder_sasa) * 100
    #else:
    #    interface_binder_fraction = 0

    # calculate surface hydrophobicity
    #binder_pose = {pose.pdb_info().chain(pose.conformation().chain_begin(i)): p for i, p in zip(range(1, pose.num_chains()+1), pose.split_by_chain())}[binder_chain]

    #layer_sel = pr.rosetta.core.select.residue_selector.LayerSelector()
    #layer_sel.set_layers(pick_core = False, pick_boundary = False, pick_surface = True)
    #surface_res = layer_sel.apply(binder_pose)

    #exp_apol_count = 0
    #total_count = 0
    #
    ## count apolar and aromatic residues at the surface
    #for i in range(1, len(surface_res) + 1):
    #    if surface_res[i] == True:
    #        res = binder_pose.residue(i)

    #        # count apolar and aromatic residues as hydrophobic
    #        if res.is_apolar() == True or res.name() == 'PHE' or res.name() == 'TRP' or res.name() == 'TYR':
    #            exp_apol_count += 1
    #        total_count += 1

    #surface_hydrophobicity = exp_apol_count/total_count

    # output interface score array and amino acid counts at the interface
    interface_scores = {
    'binder_score': 0,
    'surface_hydrophobicity': 0,
    'interface_sc': 0,
    'interface_packstat': 0,
    'interface_dG': 0,
    'interface_dSASA': 0,
    'interface_dG_SASA_ratio': 0,
    'interface_fraction': 0,
    'interface_hydrophobicity': interface_hydrophobicity,
    'interface_nres': interface_nres,
    'interface_interface_hbonds': 0,
    'interface_hbond_percentage': 0,
    'interface_delta_unsat_hbonds': 0,
    'interface_delta_unsat_hbonds_percentage':0
    }

    # round to two decimal places
    interface_scores = {k: round(v, 2) if isinstance(v, float) else v for k, v in interface_scores.items()}

    return interface_scores, interface_AA, interface_residues_pdb_ids_str

# align pdbs to have same orientation
def align_pdbs(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    ref = load_structure(reference_pdb)
    ref = ref[ref.chain_id == reference_chain_id]
    aln = load_structure(align_pdb)
    aln = aln[aln.chain_id == align_chain_id]
    alned,_,_,_ = superimpose_homologs(ref, aln)
    save_structure(align_pdb, alned)

# calculate the rmsd without alignment
def unaligned_rmsd(reference_pdb, align_pdb, reference_chain_id, align_chain_id):
    ref = load_structure(reference_pdb)
    ref = ref[ref.chain_id == reference_chain_id]
    aln = load_structure(align_pdb)
    aln = aln[aln.chain_id == align_chain_id]
    try:
        if ref.shape > aln.shape:
            return round(rmspd(aln, ref), 2)
        return round(rmspd(ref, aln), 2)
    except Exception as e:
        print(e, reference_pdb, reference_chain_id, align_pdb, align_chain_id)
        with open('errors.txt', 'a') as f:
            f.writelines('\n'.join([str(e), reference_pdb, reference_chain_id, align_pdb, align_chain_id])+'\n')
        return -1

# Relax designed structure
def pr_relax(pdb_file, relaxed_pdb_path):
    if not os.path.exists(relaxed_pdb_path):
        # Generate pose
        pose = pr.pose_from_pdb(pdb_file)
        start_pose = pose.clone()

        ### Generate movemaps
        mmf = MoveMap()
        mmf.set_chi(True) # enable sidechain movement
        mmf.set_bb(True) # enable backbone movement, can be disabled to increase speed by 30% but makes metrics look worse on average
        mmf.set_jump(False) # disable whole chain movement

        # Run FastRelax
        fastrelax = FastRelax()
        scorefxn = pr.get_fa_scorefxn()
        fastrelax.set_scorefxn(scorefxn)
        fastrelax.set_movemap(mmf) # set MoveMap
        fastrelax.max_iter(200) # default iterations is 2500
        fastrelax.min_type("lbfgs_armijo_nonmonotone")
        fastrelax.constrain_relax_to_start_coords(True)
        fastrelax.apply(pose)

        # Align relaxed structure to original trajectory
        align = AlignChainMover()
        align.source_chain(0)
        align.target_chain(0)
        align.pose(start_pose)
        align.apply(pose)

        # Copy B factors from start_pose to pose
        for resid in range(1, pose.total_residue() + 1):
            if pose.residue(resid).is_protein():
                # Get the B factor of the first heavy atom in the residue
                bfactor = start_pose.pdb_info().bfactor(resid, 1)
                for atom_id in range(1, pose.residue(resid).natoms() + 1):
                    pose.pdb_info().bfactor(resid, atom_id, bfactor)

        # output relaxed and aligned PDB
        pose.dump_pdb(relaxed_pdb_path)
        clean_pdb(relaxed_pdb_path)


MODRES = {'MSE':'MET','MLY':'LYS','FME':'MET','HYP':'PRO',
          'TPO':'THR','CSO':'CYS','SEP':'SER','M3L':'LYS',
          'HSK':'HIS','SAC':'SER','PCA':'GLU','DAL':'ALA',
          'CME':'CYS','CSD':'CYS','OCS':'CYS','DPR':'PRO',
          'B3K':'LYS','ALY':'LYS','YCM':'CYS','MLZ':'LYS',
          '4BF':'TYR','KCX':'LYS','B3E':'GLU','B3D':'ASP',
          'HZP':'PRO','CSX':'CYS','BAL':'ALA','HIC':'HIS',
          'DBZ':'ALA','DCY':'CYS','DVA':'VAL','NLE':'LEU',
          'SMC':'CYS','AGM':'ARG','B3A':'ALA','DAS':'ASP',
          'DLY':'LYS','DSN':'SER','DTH':'THR','GL3':'GLY',
          'HY3':'PRO','LLP':'LYS','MGN':'GLN','MHS':'HIS',
          'TRQ':'TRP','B3Y':'TYR','PHI':'PHE','PTR':'TYR',
          'TYS':'TYR','IAS':'ASP','GPL':'LYS','KYN':'TRP',
          'CSD':'CYS','SEC':'CYS'}

def pdb_to_string(pdb_file, chains=None, models=[1]):
  '''read pdb file and return as string'''

  if chains is not None:
    if "," in chains: chains = chains.split(",")
    if not isinstance(chains,list): chains = [chains]
  if models is not None:
    if not isinstance(models,list): models = [models]

  modres = {**MODRES}
  lines = []
  seen = []
  model = 1
  for line in open(pdb_file,"rb"):
    line = line.decode("utf-8","ignore").rstrip()
    if line[:5] == "MODEL":
      model = int(line[5:])
    if models is None or model in models:
      if line[:6] == "MODRES":
        k = line[12:15]
        v = line[24:27]
        if k not in modres and v in residue_constants.restype_3to1:
          modres[k] = v
      if line[:6] == "HETATM":
        k = line[17:20]
        if k in modres:
          line = "ATOM  "+line[6:17]+modres[k]+line[20:]
      if line[:4] == "ATOM":
        chain = line[21:22]
        if chains is None or chain in chains:
          atom = line[12:12+4].strip()
          resi = line[17:17+3]
          resn = line[22:22+5].strip()
          if resn[-1].isalpha(): # alternative atom
            resn = resn[:-1]
            line = line[:26]+" "+line[27:]
          key = f"{model}_{chain}_{resn}_{resi}_{atom}"
          if key not in seen: # skip alternative placements
            lines.append(line)
            seen.append(key)
      if line[:5] == "MODEL" or line[:3] == "TER" or line[:6] == "ENDMDL":
        lines.append(line)
  return "\n".join(lines)

def relax_protein(pdb_file: str, relaxed_pdb_path: str, max_iterations: int = 2000, tolerance: float = 2.39, stiffness: float = 1.0):
    if not os.path.exists(relaxed_pdb_path):
        pdb_str = pdb_to_string(pdb_file)
        protein_obj = protein.from_pdb_string(pdb_str)
        amber_relaxer = relax.AmberRelaxation(
            max_iterations=max_iterations,
            tolerance=tolerance,
            stiffness=stiffness,
            exclude_residues=[],
            max_outer_iterations=3,
            use_gpu=True
        )
        relaxed_pdb_lines, debug_data, _ = amber_relaxer.process(prot=protein_obj)
        with open(relaxed_pdb_path, 'w') as f:
            f.write(relaxed_pdb_lines)
        chain = protein.PDB_CHAIN_IDS[protein_obj.chain_index]
        align_pdbs(pdb_file, relaxed_pdb_path, chain, chain)
        return debug_data['final_energy']