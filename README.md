#ENCODE to Ecosystem Data Hub converter
These scripts fetch ENCODE metadata from the ENCODE portal, and output it in the Ecosystem Data Hub JSON format. Use fetch_all_exp_jsons.py to download datasets for all experiments, or individual experiment scripts. 

### Usage

#### Parameters

* assembly: Assembly name (e.g. hg19, hg38, mm10)
* taxon-id: Species taxonomy id (e.g. 9606 for human)

#### Example
```
cd IHEC_json_converter
python ./fetch_all_exp_jsons.py --assembly=hg19 --taxon-id=9606
```
