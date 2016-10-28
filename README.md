#ENCODE to Ecosystem Data Hub converter
These scripts fetch ENCODE metadata from the ENCODE portal, and output it in the Ecosystem Data Hub JSON format. Use fetch_all_exp_jsons.py to download datasets for all experiments, or individual experiment scripts. 

### Usage

#### Parameters

* assembly: Assembly name (e.g. hg19, GRCh38, mm10)
* taxon-id: Species taxonomy id (e.g. 9606 for human, 10090 for mouse)

#### Example
```
cd IHEC_json_converter
python3 ./fetch_all_exp_jsons.py --assembly=hg19 --taxon-id=9606
```

### Credits

Original development by ENCODE team, available here: https://github.com/kpaskov/PaskovEncodeScripts 