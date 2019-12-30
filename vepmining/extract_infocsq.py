import requests, sys, json ,  argparse 

freq=[]
server="https://rest.ensembl.org"
ext = "/vep/human/id/"+ "rs36080947" + "?"
r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
if not r.ok:
    r.raise_for_status()
    sys.exit()
    print(r)

decoded = r.json()
#pprint.pprint(decoded)

freq.append(decoded[0]["colocated_variants"][0]["frequencies"])
freq.append(decoded[0]["most_severe_consequence"])



