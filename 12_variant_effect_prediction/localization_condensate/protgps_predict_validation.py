import sys
import os
import pickle
from argparse import Namespace

import torch

PROTGPS_DIR = "/u/project/kappel/ddcohn/localization_tools/protgps"
sys.path.append(PROTGPS_DIR)
from protgps.utils.loading import get_object

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

COMPARTMENT_CLASSES = [
    "nuclear_speckle",
    "p-body",
    "pml-bdoy",
    "post_synaptic_density",
    "stress_granule",
    "chromosome",
    "nucleolus",
    "nuclear_pore_complex",
    "cajal_body",
    "rna_granule",
    "cell_junction",
    "transcriptional",
]


def load_model(snargs):
    modelpath = snargs.model_path
    model = get_object(snargs.lightning_name, "lightning")(snargs)
    model = model.load_from_checkpoint(
        checkpoint_path=modelpath,
        strict=not snargs.relax_checkpoint_matching,
        **{"args": snargs},
    )
    return model


@torch.no_grad()
def predict_condensates(model, sequences, batch_size=1, round_=True):
    scores = []
    for i in range(0, len(sequences), batch_size):
        batch = sequences[i:i + batch_size]
        out = model.model({"x": batch})
        s = torch.sigmoid(out["logit"]).to("cpu")
        scores.append(s)
    scores = torch.vstack(scores)
    if round_:
        scores = torch.round(scores, decimals=3)
    return scores


args = Namespace(**pickle.load(open(
    os.path.join(PROTGPS_DIR, "checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729.args"), "rb")))
args.model_path = os.path.join(
    PROTGPS_DIR, "checkpoints/protgps/32bf44b16a4e770a674896b81dfb3729epoch=26.ckpt")
args.pretrained_hub_dir = "/u/project/kappel/ddcohn/model_cache/torch_hub"
os.makedirs(args.pretrained_hub_dir, exist_ok=True)

model = load_model(args)
model.eval()
model = model.to(device)

sequences = [
    # UniProt O15116
    "MNYMPGTASLIEDIDKKHLVLLRDGRTLIGFLRSIDQFANLVLHQTVERIHVGKKYGDIPRGIFVVRGENVVLLGEIDLEKESDTPLQQVSIEEILEEQRVEQQTKLEAEKLKVQALKDRGLSIPRADTLDEY",
    # UniProt P38432
    "MAASETVRLRLQFDYPPPATPHCTAFWLLVDLNRCRVVTDLISLIRQRFGFSSGAFLGLYLEGGLLPPAESARLVRDNDCLRVKLEERGVAENSVVISNGDINLSLRKAKKRAFQLEEGEETEPDCKYSKKHWKSRENNNNNEKVLDLEPKAVTDQTVSKKNKRKNKATCGTVGDDNEEAKRKSPKKKEKCEYKKKAKNPKSPKVQAVKDWANQRCSSPKGSARNSLVKAKRKGSVSVCSKESPSSSSESESCDESISDGPSKVTLEARNSSEKLPTELSKEEPSTKNTTADKLAIKLGFSLTPSKGKTSGTTSSSSDSSAESDDQCLMSSSTPECAAGFLKTVGLFAGRGRPGPGLSSQTAGAAGWRRSGSNGGGQAPGASPSVSLPASLGRGWGREENLFSWKGAKGRGMRGRGRGRGHPVSCVVNRSTDNQRQQQLNDVVKNSSTIIQNPVETPKKDYSLLPLLAAAPQVGEKIAFKLLELTSSYSPDVSDYKEGRILSHNPETQQVDIEILSSLPALREPGKFDLVYHNENGAEVVEYAVTQESKITVFWKELIDPRLIIESPSNTSSTEPA",
]

scores = predict_condensates(model, sequences, batch_size=1)

for seq_i, seq in enumerate(sequences):
    print(f"--- sequence {seq_i} ---")
    for j, condensate in enumerate(COMPARTMENT_CLASSES):
        print(f"  {condensate.upper()}_Score {scores[seq_i, j].item()}")

print("\nExpected (from notebook reference):")
print("O15116: p-body=0.999, nucleolus=0.001, nuclear_pore_complex=0.001, rest 0.0")
print("P38432: pml-bdoy=0.003, nucleolus=0.004, nuclear_pore_complex=0.002, cajal_body=0.992, rest 0.0")
