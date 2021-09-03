import os
import sys
import time
import requests
from oss2 import headers
from oss2 import SizedFileAdapter, determine_part_size
from oss2.models import PartInfo
import json
import uuid
import argparse
import oss2
import numpy as np

def read_in_chunks(f, blocksize=1024, chunks=-1, tot=None):
    """Lazy function (generator) to read a file piece by piece.
    Default chunk size: 1k."""
    now = 0
    while chunks:
        data = f.read(blocksize)
        if not data:
            break
        yield data
        now += blocksize
        if tot is not None and now % 102400 < blocksize:
            print(("Has uploaded %.2f" % (now/tot*100)) + "%", end="\r")
        chunks -= 1


def is_float(n):
    try:
        float(n)
    except ValueError:
        return False
    return True


class Hermite():
    def __init__(self, url, email=None, password=None):
        self.url = url
        self.endpoint = "oss-cn-zhangjiakou.aliyuncs.com"
        self.bucket = "dp-tech-zhangjiakou"
        self.ossurl = "https://" + self.bucket + "." + self.endpoint + "/"
        self.sess = requests.session()
        self.login_info = None
        if email is not None:
            login_info = self.login(email, password)
            self.login_info = login_info

    def post(self, url, **kwargs):
        for i in range(5):
            try:
                print(url)
                #data = kwargs.get('json')
                # r = self.sess.post(url, json=data)
                return self.sess.post(url, **kwargs)
            except requests.exceptions.ConnectionError:
                time.sleep(1)
                continue
        raise requests.exceptions.ConnectionError

    def login(self, email, password):
        # print("Login..")
        data = {"email": email, "password": password}
        # res = self.post(self.url + "/account/login", json=data)
        res = self.post(self.url + "/Login", json=data)
        print("Respose:",res)
        print(json.loads(res.text)['code'])
        if json.loads(res.text)["code"] != 200:
            raise KeyError("Cannot login to this account!")
        json.dump(json.loads(res.text)["data"]["token"], open("temp", "w"))
        return json.loads(res.text)

    def upload_file_to_oss(self, filepath, OSSpath, bucket):
        total_size = os.path.getsize(filepath)
        if total_size < 1024**3*5: # File size is less than 5G.
            bucket.put_object_from_file(OSSpath, filepath)
        else:
            part_size = determine_part_size(total_size)
            upload_id = bucket.init_multipart_upload(OSSpath).upload_id
            parts = []
            with open(filepath, 'rb') as fileobj:
                part_number = 1
                offset = 0
                while offset < total_size:
                    num_to_upload = min(part_size, total_size - offset)
                    result = bucket.upload_part(OSSpath, upload_id, part_number,
                                                SizedFileAdapter(fileobj, num_to_upload))
                    parts.append(PartInfo(part_number, result.etag))
                    offset += num_to_upload
                    part_number += 1
            bucket.complete_multipart_upload(OSSpath, upload_id, parts)
        return OSSpath

    def upload_to_oss(self, target_file):
        file_uuid = uuid.uuid1().hex
        oss_path = "Hermite/upload/" + file_uuid + os.path.splitext(target_file)[1]
        bucket = self.get_oss_bucket()
        self.upload_file_to_oss(target_file, oss_path, bucket)
        oss_path = self.ossurl + oss_path
        print("File " + target_file + " uploaded to " + oss_path)
        return oss_path


class Hermite_func:
    def __init__(self, params, jobname="Hermite", endpoint="oss-cn-zhangjiakou.aliyuncs.com", bucket_name="dp-tech-zhangjiakou"):
        self.params = params
        self.hermite = hermite_login(params.get("pass_confirm", False))
        self.jobname = self.params.get("jobname", jobname)
        self.endpoint = self.params.get("endpoint", endpoint)
        self.bucket_name = self.params.get("bucket_name", bucket_name)
        self.token = self.hermite.login_info["data"]["token"]


class submit_file(Hermite_func):
    def __init__(self, params):
        super(submit_file, self).__init__(params, "submit_file")

    def upload_file(self):
        file = self.params.get("file", None)
        if file is None:
            raise ValueError("Please offer a protein file.")
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": os.path.split(file)[-1], "type": self.params.get("type", "")}
        )
        data = json.loads(res.text)["data"]
        print(data["oss_url"])
        file_url = data["oss_url"] + '/' + os.path.split(file)[-1]
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            file, 
            file_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        print(file_url)


class protein_preparation(Hermite_func):
    def __init__(self, params):
        super(protein_preparation, self).__init__(params, "protprep")
        self.download_path = self.params.get("download_path", None)
        if self.download_path:
            os.makedirs(self.download_path, exist_ok=True)

    def upload_protein(self):
        pdbfile = self.params.get("protein_file", None)
        if pdbfile is None:
            raise ValueError("Please offer a protein file.")
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": os.path.split(pdbfile)[-1], "type": "protein"}
        )
        print(res.text)
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + '/' + os.path.split(pdbfile)[-1]
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.url = data["oss_url"]#protein_url

    def select_group(self):
        group = self.params.get("select_group", None)
        if group is None:
            res = self.hermite.post(
                self.hermite.url + "" + "/DetectStructGroups", 
                json={"struct_oss_path": self.url}, 
                headers={'Authorization': 'jwt ' + self.token}
            )
            
            ### select group
            print("## Hermite ##\tselect groups.")
            print(res.text)
            group_info = json.loads(res.text)["data"]["group_info"]          
            print("## Hermite ##\tgroups of protein:")
            print("## Hermite ##\t|---------------------------------------------")
            print("## Hermite ##\t| idx\t| chain type\t| chain\t| start-end")
            print("## Hermite ##\t|---------------------------------------------")
            for i in group_info:
                print("## Hermite ##\t| %s"%("\t| ".join([str(p) for p in i[:]])))
            print("## Hermite ##\t|---------------------------------------------")
            group = input("## Hermite ##\tplease enter the group(s) to keep (split by ,): ")
            group = [int(i) for i in group.split(",")]
            print("## Hermite ##\tgroup(s) %s has been selected."%(", ".join([str(g) for g in group])))
        res = self.hermite.post(
            self.hermite.url + "" + "/SelectStructGroups", 
            json={"struct_oss_path": self.url, "group_select": group}, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)
        self.url = json.loads(res.text)["data"]["url_select_protein"]
        self.herero_url = json.loads(res.text)["data"]["url_select_hetero"]

    def repair(self):
        try:
            list_missing = self.params.get("list_missing", '').split(",")
        except:
            list_missing = None
        if not list_missing:
            res = self.hermite.post(
                self.hermite.url + "" + "/DetectProteinMissing", 
                json={"protein_oss_path": self.url}, 
                headers={'Authorization': 'jwt ' + self.token}
            )
            missing_AA = json.loads(res.text)["data"]["missing_residue_info"]
            if missing_AA == []:
                print("## Hermite ##\tthere is no missing part.")
            else:
                print("## Hermite ##\tthere are %d missing parts."%len(missing_AA))
                print("## Hermite ##\t|---------------------------------------------")
                print("## Hermite ##\t| idx\t| chain\t| s-e\t| sequence")
                print("## Hermite ##\t|---------------------------------------------")
                for midx, mc, mse, ms in missing_AA:
                    print("## Hermite ##\t|%d\t| %s\t| %s\t| %s"%(midx, mc, mse, ", ".join(ms)))
                print("## Hermite ##\t|---------------------------------------------")
                list_missing = input("## Hermite ##\tplease enter the missing part(s) to repair (split by ,): ")
                list_missing = list_missing.split(",")
        res = self.hermite.post(
            # self.hermite.url + "" + "/RepairProtein", 
            self.hermite.url + "" + "/RepairAndPrepareProtein", 
            json={
                "protein_oss_path": self.url, 
                "hereo_oss_path" :self.herero_url, 
                "missing_repair": list_missing,
                "add_missing_side_chains": self.params["add_missing_side_chains"],
                "add_hydrogens": self.params["add_hydrogens"],
                "optimize_the_hydrogen_bonding_network": self.params["optimize_the_hydrogen_bonding_network"],
                "protonation_state": self.params["protonation_state"],
                "protonation_ph": self.params["protonation_ph"],
            }, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        self.url = json.loads(res.text)["data"]["url"]
        if self.download_path:
            _download_path = os.path.join(self.download_path, "protein.pdb")
            print("## Hermite ##\tdownload prepared protein(pdb) to %s" % _download_path)
            os.system("wget %s -O %s" % (self.url, _download_path))
        print(json.loads(res.text))

    def protein_prep(self):
        print("No adt process:",self.url)
        print("\n\nPortein_prep\n\n")
        res = self.hermite.post(
            self.hermite.url + "" + "/PrepareProtein", 
            json={"protein_oss_path": self.url}, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print("\n\nPortein_prep\n\n")
        print(res.text) 
        self.url = json.loads(res.text)["data"]["protein_url"]
        print("## Hermite ##\tyour prepared protein oss path is:")
        print("## Hermite ##\t%s"%self.url)
        print("## Hermite ##\tplease save the url which will be used when submit docking job.")
        
        if self.download_path:
            _download_path = os.path.join(self.download_path, "protein.pdbqt")
            print("## Hermite ##\tdownload prepared protein(pdbqt) to %s"%_download_path)
            os.system("wget %s -O %s"%(self.url, _download_path))

    def set_flexible(self):
        try:
            flex_res = self.params.get("flex_res", "").split(",")
        except:
            flex_res = []
        
        center_x = self.params.get("center_x", None)
        center_y = self.params.get("center_y", None)
        center_z = self.params.get("center_z", None)
        cutoff = self.params.get("cutoff", None)
        res = self.hermite.post(
            self.hermite.url + "" + "/GetCenterData", 
            json={"protein_oss_path": self.url}, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        center_dict = json.loads(res.text)["data"]["res_info"]    
        center_data = [[
            cd["center_x"],
            cd["center_y"],
            cd["center_z"],
            cd["residue_name"],
            idx
            ] for idx, cd in enumerate(center_dict)]
        residue_dict = {int(cd["residue_name"][3:]):idx for idx, cd in enumerate(center_dict)}
        v_residue_dict = {idx:int(cd["residue_name"][3:]) for idx, cd in enumerate(center_dict)}
        residue_name_dict = {idx:cd["residue_name"][:3] for idx, cd in enumerate(center_dict)}
        if flex_res == []:
            print("## Hermite ##\tthere are %d residues in the protein:"%len(center_data))
            print("## Hermite ##\t|---------------------------------------------")
            print("## Hermite ##\t| idx\t| res\t| x\t| y\t| z")
            print("## Hermite ##\t|---------------------------------------------")
            for i in center_data:
                print("## Hermite ##\t| %s\t| %s\t| %.1f\t| %.1f\t| %.1f"%(i[3][3:],i[3][:3],i[0],i[1],i[2]))
            print("## Hermite ##\t|---------------------------------------------")
            flex_res_idx = []
            if self.params.get("distance_select", False):
                print("## Hermite ##\tselect flexible residues based on distance.")
                if center_x is None:
                    center_x = float(input("## Hermite ##\tplease enter the coordinate x of box center: "))
                else:
                    print("## Hermite ##\tpreset coordinate x of box center: %.1f"%center_x)
                if center_y is None:
                    center_y = float(input("## Hermite ##\tplease enter the coordinate y of box center: "))
                else:
                    print("## Hermite ##\tpreset coordinate y of box center: %.1f"%center_y)
                if center_z is None:
                    center_z = float(input("## Hermite ##\tplease enter the coordinate z of box center: "))
                else:
                    print("## Hermite ##\tpreset coordinate z of box center: %.1f"%center_z)
                cd = np.array(center_data)[:, :3].astype(float)
                distance = cd - np.array([center_x,center_y,center_z]).reshape(1,3)
                distance = np.sqrt((distance**2).sum(-1))
                cutoff = self.params.get("cutoff", None)
                if cutoff is None:
                    cutoff = float(input("## Hermite ##\tplease enter the cutoff: "))
                else:
                    print("## Hermite ##\tpreset cutoff: %.1f"%cutoff)
                while True:
                    mask = distance < cutoff
                    flex_res_idx = np.array(center_data)[mask][:,-1].astype(int).tolist()
                    flex_res_name = np.array(center_data)[mask][:,-2].tolist()
                    flex_res_dist = distance[mask].tolist()
                    print("## Hermite ##\tresidues within cutoff %.1f (%d): "%(cutoff, len(flex_res_idx)))
                    print("## Hermite ##\t|---------------------------------------------")
                    print("## Hermite ##\t| idx\t| res\t| distance")
                    print("## Hermite ##\t|---------------------------------------------")
                    temp_idx_name_dist = [[n[3:],n[:3],d] for n,d in zip(flex_res_name, flex_res_dist)]
                    temp_idx_name_dist.sort(key=lambda x:x[-1])
                    for n1, n2, d in temp_idx_name_dist:
                        print("## Hermite ##\t| %s\t| %s\t| %.1f"%(n1,n2,d))
                    print("## Hermite ##\t|---------------------------------------------")
                    if input("## Hermite ##\treset cutoff? (Y/N) ").lower() == "n":
                        break
                    cutoff = float(input("## Hermite ##\tplease enter the cutoff: "))
                mask = distance < cutoff
                flex_res_idx = np.array(center_data)[mask][:,-1].astype(int).tolist()
                flex_res_name = np.array(center_data)[mask][:,-2].tolist()
                flex_res_dist = distance[mask].tolist()
                flex_res_idx = np.array(center_data)[mask][:,-1].astype(int).tolist()
                print("## Hermite ##\tresidues within cutoff %.1f (%d): "%(cutoff, len(flex_res_idx)))
                print("## Hermite ##\t|---------------------------------------------")
                print("## Hermite ##\t| idx\t| res\t| distance")
                print("## Hermite ##\t|---------------------------------------------")
                temp_idx_name_dist = [[n[3:],n[:3],d] for n,d in zip(flex_res_name, flex_res_dist)]
                temp_idx_name_dist.sort(key=lambda x:x[-1])
                for n1, n2, d in temp_idx_name_dist:
                    print("## Hermite ##\t| %s\t| %s\t| %.1f"%(n1,n2,d))
                print("## Hermite ##\t|---------------------------------------------")
            flex_res_idx.extend([residue_dict[int(ni)] for ni in flex_res])
            flex_res = []
            temp_flex_res_idx = flex_res_idx + []
            while True:
                temp_flex_res_idx.sort()
                print("## Hermite ##\tcurrent flexible residues (%d):"%len(temp_flex_res_idx))
                print("## Hermite ##\t|---------------------------------------------")
                print("## Hermite ##\t| idx\t| res")
                print("## Hermite ##\t|---------------------------------------------")
                for idx, res in zip(
                    [v_residue_dict[p] for p in temp_flex_res_idx],
                    [residue_name_dict[p] for p in temp_flex_res_idx]
                    ):
                    print("## Hermite ##\t| %d\t| %s"%(idx, res))
                print("## Hermite ##\t|---------------------------------------------")
                if input("## Hermite ##\tcontinue to add? (Y/N) ").lower() == "n":
                    break
                temp_flex_res = ""
                temp_flex_res = input("## Hermite ##\tplease enter the index of residues to be flexible (additionally): ")
                temp_flex_res = temp_flex_res.split(",")
                flex_res.extend(temp_flex_res)
                temp_flex_res_idx.extend([residue_dict[int(ni)] for ni in temp_flex_res])
        try:
            flex_res_idx.extend([residue_dict[int(ni)] for ni in flex_res])
        except:
            flex_res_idx = []
            flex_res_idx.extend([residue_dict[int(ni)] for ni in flex_res])
        print("## Hermite ##\t| final flexible residues (%d):"%len(flex_res_idx))
        print("## Hermite ##\t|---------------------------------------------")
        print("## Hermite ##\t| idx\t| res")
        print("## Hermite ##\t|---------------------------------------------")
        temp_idx_name = [[v_residue_dict[p], residue_name_dict[p]] for p in flex_res_idx]
        temp_idx_name.sort(key=lambda x:x[0])
        for idx, res in temp_idx_name:
            print("## Hermite ##\t| %d\t| %s"%(idx, res))
        print("## Hermite ##\t|---------------------------------------------")
        res = self.hermite.post(
            self.hermite.url + "" + "/FlexPrep", 
            json={"protein_oss_path": self.url, "flex_res_idx": flex_res_idx}, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        
        print("## Hermite ##\tyour prepared protein oss path is:")
        print("## Hermite ##\t%s"%self.url)
        print("## Hermite ##\tplease save the url which will be used when submit docking job.")
        if self.download_path:
            protein_flex_url = json.loads(res.text)["data"]["protein_flex_url"]
            protein_rigid_url = json.loads(res.text)["data"]["protein_rigid_url"]
            _download_path = os.path.join(self.download_path, "protein_flex.pdbqt")
            print("## Hermite ##\tdownload prepared flex protein(pdbqt) to %s"%_download_path)
            os.system("wget %s -O %s"%(protein_flex_url, _download_path))
            _download_path = os.path.join(self.download_path, "protein_rigid.pdbqt")
            print("## Hermite ##\tdownload prepared rigid protein(pdbqt) to %s"%_download_path)
            os.system("wget %s -O %s"%(protein_rigid_url, _download_path))


class MD(Hermite_func):
    def __init__(self, params):
        super(MD, self).__init__(params, "MD")

    def upload_protein(self):
        try:
            complex_pdb = self.params.get("complex_pdb", None).split(",")
            del self.params["complex_pdb"]
        except:
            complex_pdb = None
        if complex_pdb is None:
            raise ValueError("Please offer at least one complex file.")
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": self.jobname, "type": "protein"}
            )
        data = json.loads(res.text)["data"]
        oss_url = data["oss_url"]
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading complex file...")
        for cpdb in complex_pdb:
            protein_url = os.path.join(oss_url, os.path.split(cpdb)[-1])
            self.hermite.upload_file_to_oss(
                cpdb, 
                protein_url[len(self.hermite.ossurl):], 
                bucket
                )
        self.url = oss_url
        print("## Hermite ##\tuploading sucessed.")

    def submitMD(self):
        dicts = {q:self.params[q] for q in self.params.keys()}
        dicts["complex_oss_path"] = self.url
        res = self.hermite.post(
            self.hermite.url + "" + "/SubmitMD", json=dicts, 
            headers={'Authorization': 'jwt ' + self.token}
            )
        print(res.text)


class ligand_preparation(Hermite_func):
    def __init__(self, params):
        super(ligand_preparation, self).__init__(params, "ligprep")

    def upload_ligand(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": self.jobname, "type": "ligand"}
        )
        data = json.loads(res.text)["data"]
        self.params["ligand_url"] = data["oss_url"]
        ligand_files = self.params.get("ligfiles", "").split(",")
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        for lf in ligand_files:
            lu = self.params["ligand_url"][len(self.hermite.ossurl):] + "/" + os.path.split(lf)[-1]
            self.hermite.upload_file_to_oss(lf, lu, bucket)

    def ligand_prep(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/SubmitLigandPrep", json=self.params, 
            headers={'Authorization': 'jwt ' + self.token}
            )
        print(res.text)


class submit_dock(Hermite_func):
    def __init__(self, params):
        super(submit_dock, self).__init__(params, "dock")

    def upload_config(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/UploadDockConfig",
            json={
                "name": self.params.get("jobname", "dock_config"),
                "center_x": self.params["center_x"],
                "center_y": self.params["center_y"],
                "center_z": self.params["center_z"],
                "size_x": self.params["size_x"],
                "size_y": self.params["size_y"],
                "size_z": self.params["size_z"],
                "spacing": self.params["spacing"],
            }, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)
        self.params["config_oss_path"] = json.loads(res.text)["data"]["config_url"]

    def upload_protein(self):
        pdbfile = self.params.get("protein_file", None)
        if pdbfile is None:
            raise ValueError("Please offer a protein file.")
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": os.path.split(pdbfile)[-1], "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + '/' + os.path.split(pdbfile)[-1]
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.params["protein_oss_path"] = data["oss_url"]

    def upload_ligand(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": self.jobname, "type": "ligand"}
        )
        data = json.loads(res.text)["data"]
        self.params["ligand_oss_path"] = data["oss_url"]
        ligand_file = self.params.get("ligfile", "")
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        lu = self.params["ligands_oss_path"][len(self.hermite.ossurl):] + "/" + os.path.split(ligand_file)[-1]
        self.hermite.upload_file_to_oss(ligand_file, lu, bucket)

    def submit_dock(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/SubmitDock",
            json=self.params, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)


class ligand_binding_site_prediction(Hermite_func):
    def __init__(self, params):
        super(ligand_binding_site_prediction, self).__init__(params, "ligand_binding_site_prediction")
    
    def upload_protein(self):
        pdbfile = self.params.get("protein_file", None)
        if pdbfile is None:
            raise ValueError("Please offer a protein file.")
        pdbname = os.path.split(pdbfile)[-1].split(".")[0]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": pdbname + ".pdb", "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + "/" + pdbname + ".pdb"
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.url = data["oss_url"]#protein_url

    def submitLBS(self):
        res = self.hermite.post(
            self.hermite.url + "" + "/SubmitLBS", 
            json={"protein_oss_path": self.url}, 
            headers={'Authorization': 'jwt ' + self.token}
            )
        print(res.text)

class ADMET(Hermite_func):
    def __init__(self, params):
        super(ADMET, self).__init__(params, "ADMET")
    
    def upload_files(self):
        file_list = self.params.get("files", None)
        if not file_list:
            raise ValueError("No input file.")
        file_list = file_list.split(',')
        self.url_list = list()
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": self.jobname, "type": "ligand"}
        )
        data = json.loads(res.text)["data"]
        oss_url = data["oss_url"]
        auth = oss2.StsAuth(data["AccessKeyId"],
                    data["AccessKeySecret"],
                    data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        for file in file_list:
            filename = os.path.split(file)[-1]
            file_path = os.path.join(oss_url[len(self.hermite.ossurl):], filename)
            print("## Hermite ##\tuploading ligand file...")
            self.hermite.upload_file_to_oss(
                file, 
                file_path, 
                bucket
            )
            self.url_list.append(os.path.join(oss_url, filename))
        
    def submitADMET(self):
        dicts = {q: self.params[q] for q in self.params.keys()}
        dicts["oss_path_list"] = self.url_list
        outtypes = self.params.get("output_types", None)
        outtypes = outtypes.split(',')
        dicts['props'] = outtypes
        res = self.hermite.post(
            self.hermite.url + "" + "/SubmitADMET", 
            json=dicts, 
            headers={'Authorization': 'jwt ' + self.token}
            )
        print(res.text)

class LoopModelling(Hermite_func):
    def __init__(self, params):
        super(LoopModelling, self).__init__(params, "LoopModelling")
        self.seq_url = None
        self.protein_url = None
    
    def upload_file(self):
        pdbfile = self.params.get("pdbfile", None)
        if pdbfile is None:
            raise FileNotFoundError("Please offer a protein file.")
        seqfile = self.params.get("seqfile", None)
        if seqfile is None:
            print("Could not find sequence file, protein file must have sequence information!")
        
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": "loop_modelling", "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        pdbname = os.path.split(pdbfile)[-1]
        pdb_url = data["oss_url"] + '/' + pdbname
        auth = oss2.StsAuth(data["AccessKeyId"],
                    data["AccessKeySecret"],
                    data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading pdb file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            pdb_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        if seqfile is not None:
            seq_name = os.path.split(seqfile)[-1]
            seq_url = data["oss_url"] + '/' + seq_name
            print("## Hermite ##\tuploading seq file...")
            self.hermite.upload_file_to_oss(
                seqfile, 
                seq_url[len(self.hermite.ossurl):], 
                bucket
            )
            print("## Hermite ##\tuploading sucessed.")
        else:
            seq_url = None
        self.seq_url = seq_url
        self.protein_url = pdb_url

    def submitLoopModelling(self):
        self.params["protein_oss_path"] = self.protein_url
        self.params["seq_oss_path"] = self.seq_url
        res = self.hermite.post(self.hermite.url + "" + "/SubmitLoopModelling", 
        json=self.params, 
        headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)

class LoopOpt(Hermite_func):
    def __init__(self, params):
        super(LoopOpt, self).__init__(params, "LoopOpt")
    
    def upload_protein(self):
        pdbfile = self.params.get("pdbfile", None)
        if pdbfile is None:
            raise FileNotFoundError("Please offer a protein file.")
        pdbname = os.path.split(pdbfile)[-1]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": pdbname, "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + "/" + pdbname
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.protein_url = data["oss_url"]
    
    def submitLoopOpt(self):
        self.params["protein_oss_path"] = self.protein_url
        res = self.hermite.post(self.hermite.url + "" + "/SubmitLoopOpt", 
        json=self.params, 
        headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)
    

class StructureRefinement(Hermite_func):
    def __init__(self, params):
        super(StructureRefinement, self).__init__(params, "StructureRefinement")
    
    def upload_protein(self):
        pdbfile = self.params.get("pdbfile", None)
        if pdbfile is None:
            raise FileNotFoundError("Please offer a protein file.")
        pdbname = os.path.split(pdbfile)[-1]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": pdbname, "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + "/" + pdbname
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading protein file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.protein_url = data["oss_url"]#protein_url
    
    def submitStructureRefinement(self):
        self.params["selected_residue_index"] = [int(rind) for rind in self.params["selected_residue_index"].split(',')]
        self.params["protein_oss_path"] = self.protein_url
        res = self.hermite.post(self.hermite.url + "" + "/SubmitStructureRefinement", 
            json=self.params, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)


class HomologyModelling(Hermite_func):
    def __init__(self, params):
        super(HomologyModelling, self).__init__(params, "HomologyModelling")
        self.template_url = None
    
    def upload_sequence(self):
        seqfile = self.params.get("seqfile", None)
        if seqfile is None:
            raise FileNotFoundError("Please offer a sequence file.")
        seqname = os.path.split(seqfile)[-1]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": seqname, "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        sequence_url = data["oss_url"] + "/" + seqname
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading sequence file...")
        self.hermite.upload_file_to_oss(
            seqfile, 
            sequence_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.seq_url = data["oss_url"]

    def upload_template(self):
        pdbfile = self.params.get("tmpfile", None)
        pdbname = os.path.split(pdbfile)[-1]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": pdbname, "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + "/" + pdbname
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        print("## Hermite ##\tuploading template file...")
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        print("## Hermite ##\tuploading sucessed.")
        self.template_url = data["oss_url"]
    
    def submitHomologyModelling(self):
        self.params["sequence_oss_path"] = self.seq_url
        self.params["template_oss_path"] = self.template_url
        res = self.hermite.post(self.hermite.url + "" + "/SubmitHomologyModelling", 
            json=self.params, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)


class SimilaritySearch(Hermite_func):
    def __init__(self, params):
        super(SimilaritySearch, self).__init__(params, "SimilaritySearch")

    def uploadSimSearch(self):
        res = self.hermite.post(self.hermite.url + "" + "/SubmitSimiSearch", 
            json=self.params, 
            headers={'Authorization': 'jwt ' + self.token}
            )
        print(res.text)


class VSW(Hermite_func):
    def __init__(self, params):
        super(VSW, self).__init__(params, "VSW")
        self.task_list = []
        self.workflow = []
        self.config_oss_path = None
        self.protein_oss_path = None
        self.lig_oss_path = None

    def selectTasks(self):
        selected_list = ["ligprep", "dock", "MD"]
        print("## Hermite ##\ttasks:")
        print("## Hermite ##\t|---------------------------------------------")
        for i, t in enumerate(selected_list):
            print("## Hermite ##\t| %s: %s" % (str(i + 1), t))
        print("## Hermite ##\t|---------------------------------------------")
        tasks = input("## Hermite ##\tplease enter the tasks no. by order (split by ,): ")
        tasks = [selected_list[int(i) - 1] for i in tasks.split(",")]
        self.task_list = tasks
        print("## Hermite ##\ttasks %s has been selected."%(", ".join([t for t in tasks])))

    def upload_protein(self, pdbfile):
        if pdbfile is None:
            raise ValueError("Please offer a protein file.")
        pdbname = os.path.split(pdbfile)[-1].split(".")[0]
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": pdbname + ".pdb", "type": "protein"}
        )
        data = json.loads(res.text)["data"]
        protein_url = data["oss_url"] + "/" + pdbname + ".pdb"
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        self.hermite.upload_file_to_oss(
            pdbfile, 
            protein_url[len(self.hermite.ossurl):], 
            bucket
        )
        self.protein_oss_path = data["oss_url"]

    def upload_ligand(self, ligfiles):
        res = self.hermite.post(
            self.hermite.url + "" + "/GetOSSToken", 
            json={"name": self.jobname, "type": "ligand"}
        )
        data = json.loads(res.text)["data"]
        oss_url = data["oss_url"]
        ligand_files = ligfiles.split(",")
        auth = oss2.StsAuth(data["AccessKeyId"],
            data["AccessKeySecret"],
            data["SecurityToken"])
        bucket = oss2.Bucket(auth, self.endpoint, self.bucket_name)
        for lf in ligand_files:
            lu = oss_url[len(self.hermite.ossurl):] + "/" + os.path.split(lf)[-1]
            self.hermite.upload_file_to_oss(lf, lu, bucket)
        self.lig_oss_path = oss_url

    def upload_dock_config(self, params):
        res = self.hermite.post(
            self.hermite.url + "" + "/UploadDockConfig",
            json={
                "name": params.get("jobname", "dock_config"),
                "center_x": params["center_x"],
                "center_y": params["center_y"],
                "center_z": params["center_z"],
                "size_x": params["size_x"],
                "size_y": params["size_y"],
                "size_z": params["size_z"],
                "spacing": params["spacing"],
            }, 
            headers={'Authorization': 'jwt ' + self.token}
        )
        config_oss_path = json.loads(res.text)["data"]["config_url"]
        self.config_oss_path = config_oss_path

    def uploadParams(self):
        for i, task in enumerate(self.task_list):
            print("Please input parameters for %s" % task)
            # if task == "proprep":
            #     params = {"select_group": None, "distance_select": False, "centerx": None, 
            #         "centery": None, "centerz": None, "cutoff": None, "flex_res": None, 
            #         "download_path": "", "list_missing": None, 
            #         "not_add_missing_side_chains": True, "not_add_hydrogens": True, 
            #         "not_optimize_the_hydrogen_bonding_network": True, 
            #         "not_set_protonation_state": True, 
            #         "protonation_ph": 7, "set_flexible": False, "user_confirm": False
            #     }
            if task == "ligprep":
                if i + 1 < len(self.task_list) - 1 and self.task_list[i + 1] != "dock":
                    raise KeyError("Unsupported job order")
                if i > 0 and self.task_list[i - 1] == "dock":
                    self.lig_oss_path = ""
                else:
                    ligfiles = input("Please input ligand file paths (use ',' to seperate): \n")
                    self.upload_ligand(ligfiles)
                params = {
                    "ligand_url": self.lig_oss_path, "directly_convert": False, 
                    "not_calc_prop": True, "Lipinski": False, "not_add_hydrogen": True, 
                    "ph_mean": 7.4, "ph_std": 1.0, "not_generate_tautomers": True, 
                    "number_of_tautomers": 1, "to_3d": False, "optimize_structures": False, 
                    "not_flexible": True, "ntasks": 1, "jobname": "ligprep", "user_confirm": False
                }
                inp = input("Please input the parameters, format should be like 'param1=a,param2=b', parameters you do not input will use default value.\n Parameters are: %s: \n" % params.keys())
                inp = inp.split(",")
                for pair in inp:
                    if pair:
                        p, v = pair.split('=')
                        if p not in params:
                            print("Invalid paramaters")
                            continue
                        if v == "False":
                            v = False
                        elif v == "True":
                            v = True
                        elif v == "None":
                            v = None
                        elif v.isalnum():
                            v = int(v)
                        elif is_float(v):
                            v = float(v)
                        params[p] = v
                params["jobtype"] = "ligprep"
                self.workflow.append(params)

            elif task == "dock":
                if self.config_oss_path is None:
                    config_inp = dict()
                    req_param = input("Please input center-x, center-y, center-z, size-x, size-y, size-z(use ',' to seperate): \n")
                    req_param = [float(n) for n in req_param.split(',')]
                    if not req_param or len(req_param) != 6:
                        raise KeyError("You must input these 6 parameters")
                    spac = input("You can input spacing parameters, or use default value: \n")
                    if not spac:
                        spac = 0.375
                    else:
                        spac = float(spac)
                    req_param.append(spac)
                    for i, p in enumerate(['center_x', 'center_y', 'center_z', 'size_x', 'size_y', 'size_z', 'spacing']):
                        config_inp[p] = req_param[i]
                    self.config_oss_path = self.upload_dock_config(config_inp)
                if i == 0:
                    ligfiles = input("Please input ligand file paths (use ',' to seperate)")
                    self.lig_oss_path = self.upload_ligand("ligprep", ligfiles)
                params = {
                    "protein_oss_path": self.protein_oss_path, 
                    "ligands_oss_path": self.lig_oss_path, "engine": "vina", "ad4_gpf": None, "ad4_dpf": None, 
                    "exhausiveness": 1, "number_of_modes": 8, "energy_range": 3, 
                    "number_of_ga": 10, "maximum_number_of_evals": 2500000, 
                    "maximum_number_of_generations": 27000, "keep_number": None, 
                    "not_save_ligs_only": True, "not_calc_prop": True, "ntasks": 5, 
                    "jobname": "dock", "user_confirm": False
                }
                params["jobtype"] = "dock"
                self.workflow.append(params)

            elif task == "MD":
                if i != len(self.task_list) - 1:
                    raise KeyError("Unsupported job order")
                params = {
                    "complex_oss_path": "", "timestep": 0.002, "steps": "5000000", 
                    "protein_forcefield": "amber99sb", "ligand_forcefield": "gaff", 
                    "box_type": "triclinic", "boundary": "1.0", "heavyh": False, 
                    "ion_conc": "0.15", "threshold": None, "jobname": "MD", 
                    "user_confirm": False, "use_gpu": False, "from_dock": False
                }
                params["jobtype"] = "MD"
                self.workflow.append(params)

    def submitVSW(self):
        print(self.workflow)
        res = self.hermite.post(self.hermite.url + "" + "/SubmitVSW", 
        json={"workflow": self.workflow}, 
        headers={'Authorization': 'jwt ' + self.token}
        )
        print(res.text)


class QUERY(Hermite_func):
    def __init__(self, params):
        super(QUERY, self).__init__(params, "QUERY")
        self.id = self.params.get("id", None)
    
    def query(self):
        if self.params["type"] == "USER_JOBS":
            res = self.hermite.post(self.hermite.url + "" + "/GetUserJobs", 
                headers={'Authorization': 'jwt ' + self.token}
                )
            print(res.text) 
        elif self.params["type"] == "USER_FILES":
            res = self.hermite.post(self.hermite.url + "" + "/GetUserFiles", 
                headers={'Authorization': 'jwt ' + self.token}
                )
            print(res.text)    
        elif self.params["type"] == "RES_OSS":
            if self.id is None:
                raise ValueError("please provide job id")
            res = self.hermite.post(self.hermite.url + "" + "/GetResultsOSSPath", 
                json={"id": self.id}, 
                headers={'Authorization': 'jwt ' + self.token}
                )
            print(res.text)
        elif self.params["type"] == "PARAM":
            if self.id is None:
                raise ValueError("please provide job id")
            res = self.hermite.post(self.hermite.url + "" + "/GetInputParams", 
                json={"id": self.id}, 
                headers={'Authorization': 'jwt ' + self.token}
                )
            print(res.text)
        else:
            raise TypeError("Incorrect query type!")


def hermite_login(pass_confirm):
    sdk_path = os.path.split(os.path.realpath(__file__))[0]
    config_path = os.path.join(sdk_path, "hermite_sdk_config.json")
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
    except:
        raise ValueError("## Hermite ##\tplease login first.")
    if not pass_confirm:
        cfg = input("## Hermite ##\tlogin with email [%s]?\t(Y/N)\t"%config.get("email", None))
        for i in range(3):
            if cfg.lower() == "y":
                break
            elif cfg.lower() == "n":
                print("## Hermite ##\tplease login with your email first.")
                sys.exit()
            else:
                if i == 2:
                    print("## Hermite ##\trequest failed three times, please restart.")
                    sys.exit()
                cfg = input("## Hermite ##\tplease choose Y/N: ")
    print("## Hermite ##\tlogining...")
    hermite = Hermite(
        # url="https:/.dp.tech", 
        url="http://39.99.238.157:5011", 
        email=config.get("email", None),
        password=config.get("password", None)
        )
    print("## Hermite ##\tlogin sucessed.")
    print("## Hermite ##\twelcome [%s]!"%config.get("email", None))
    return hermite


def createProject(project_name, group_name, pass_confirm=False):
    hermite = hermite_login(pass_confirm)
    token = hermite.login_info["data"]["token"]
    res = hermite.post(hermite.url + "" + "/NewProject", 
        json={"name": project_name, "group": group_name}, 
        headers={'Authorization': 'jwt ' + token})
    print(res.text)

def rerunTask(id, pass_confirm=False):
    hermite = hermite_login(pass_confirm)
    token = hermite.login_info["data"]["token"]
    res = hermite.post(hermite.url + "" + "/RerunTask", 
        json={"id": id}, 
        headers={'Authorization': 'jwt ' + token})
    print(res.text)


def main():
    parser = argparse.ArgumentParser(prog='python Hermite_SDK.py',formatter_class=argparse.RawTextHelpFormatter)

    subparsers = parser.add_subparsers(dest="cmd", help='sub-command help')

    parser_login = subparsers.add_parser("login", help="log in")
    parser_login.add_argument("-e", "--email", type=str, required=True,
        help="Email to login.")
    parser_login.add_argument("-p", "--password", type=str, required=True,
        help="Password corresponding to the user name")
    parser_login.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_project = subparsers.add_parser("project", help="create group")
    parser_project.add_argument("-p", "--project", type=str, required=True, 
        help="group name")
    parser_project.add_argument("-g", "--group", type=str, default=None, 
        help="group name")
    parser_project.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    
    parser_rerun = subparsers.add_parser("rerun", help="create group")
    parser_rerun.add_argument("--id", type=int, required=True, 
        help="job id")
    parser_rerun.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_submit = subparsers.add_parser("submit", help="submit a file")
    parser_submit.add_argument("-f", "--file", type=str, required=True, 
        help="local file path")
    parser_submit.add_argument("-d", "--type", type=str, default='protein', 
        choices=['protein', 'ligand'], help="file type")
    parser_submit.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_proprep = subparsers.add_parser("proprep", help="protein preparation")
    parser_proprep.add_argument("-f", "--protein_file", type=str, required=True,
        help="Protein file.")
    parser_proprep.add_argument("-g", "--select_group", type=str, default=None,
        help="Group to keep. If you want to choose several groups, split by ','.")
    parser_proprep.add_argument("-d", "--distance_select", action="store_true",
        help="Calculate distance between residues and box center to determine the flexible residues.")
    parser_proprep.add_argument("-cx", "--center_x", type=float, default=None,
        help="Coordinate x of box center.")
    parser_proprep.add_argument("-cy", "--center_y", type=float, default=None,
        help="Coordinate y of box center.")
    parser_proprep.add_argument("-cz", "--center_z", type=float, default=None,
        help="Coordinate z of box center.")
    parser_proprep.add_argument("-c", "--cutoff", type=float, default=None,
        help="cutoff of distance between residues and box center.")
    parser_proprep.add_argument("-fr", "--flex_res", type=str, default=None,
        help="Index of residues to set flexible. If you want to convert several indexes, split by ','.")
    parser_proprep.add_argument("-p", "--download_path", type=str, default="",
        help="Path to save protein preparation results.")
    parser_proprep.add_argument("-lm", "--list_missing", default=None, type=str,
        help="Missing part to add. (default all)") # A sequence of fasta
    parser_proprep.add_argument("--not_add_missing_side_chains", action="store_false",
        dest="add_missing_side_chains",
        help="Add missing side chains.")
    parser_proprep.add_argument("--not_add_hydrogens", action="store_false",
        dest="add_hydrogens",
        help="Add hydrogens.")
    parser_proprep.add_argument("--not_optimize_the_hydrogen_bonding_network", action="store_false",
        dest="optimize_the_hydrogen_bonding_network",
        help="Optimize the hydrogen bonding network.")
    parser_proprep.add_argument("--not_set_protonation_state", action="store_false",
        dest="protonation_state",
        help="Set protonation state.")
    parser_proprep.add_argument("-ph", "--protonation_ph", default=7,
        help="PH of protonation state.")
    parser_proprep.add_argument("--set_flexible", action="store_true",
        help="set flexible residues.")
    parser_proprep.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    parser_proprep.add_argument("--projectid", default=None, 
        help="project ID in table")

    parser_ligprep = subparsers.add_parser("ligprep", help="ligand preparation")
    parser_ligprep.add_argument("-l", "--ligfiles", required=True, type=str,
        help="Ligand files to do ligand preparation. If you want to convert several ligand files at the same time, \
              please seperate them with ','.")
    parser_ligprep.add_argument("-d", "--directly_convert", action="store_true",
        help="Directly convert format.")
    parser_ligprep.add_argument("-noP", "--not_calc_prop", dest="calc_prop", action="store_false",
        help="Not calculate molecular properties.")
    parser_ligprep.add_argument("-L", "--Lipinski", action="store_true",
        help="Filter ligands by Lipinski's ruls.")
    parser_ligprep.add_argument("-noH", "--not_add_hydrogen", dest="add_hydrogen", action="store_false",
        help="Not do inoization.")
    parser_ligprep.add_argument("-phm", "--ph_mean", type=float, default=7.4,
        help="(min pH + max pH) / 2")
    parser_ligprep.add_argument("-phs", "--ph_std", type=float, default=1.0,
        help="(max pH - min pH) / 2")
    parser_ligprep.add_argument("-noT", "--not_generate_tautomers", dest="generate_tautomers", action="store_false",
        help="Not generate tautomers.")
    parser_ligprep.add_argument("-n", "--number_of_conformation", type=int, default=4,
        help="The number of tautomers to be generated.")
    parser_ligprep.add_argument("-n3D", "--not_to_3d", dest="to_3d", action="store_false",
        help="Not convert from 2D to 3D.")
    parser_ligprep.add_argument("-o", "--optimize_structures", action="store_true",
        help="Optimize structure.")
    parser_ligprep.add_argument("-noF", "--not_flexible", dest="flexible", action="store_false",
        help="Not set ligand flexible.")
    parser_ligprep.add_argument("-j", "--jobname", type=str, default="ligprep",
        help="The name of this job.")
    parser_ligprep.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    parser_ligprep.add_argument("--projectid", default=None, 
        help="project ID in table")

    parser_dock = subparsers.add_parser("dock", help="submit docking job")
    parser_dock.add_argument("-p", "--protein_oss_path", type=str, required=True,
        help="Prepared protein file in PDBQT format.")
    parser_dock.add_argument("-l", "--ligands_oss_path", type=str, required=True,
        help="OSS path of the zipped file of the prepared ligand file in PDBQT format.")
    parser_dock.add_argument("-e", "--engine", type=str, default="vina",
        choices=["vina", "smina", "ad4"],
        help="The engine used to dock.")
    parser_dock.add_argument("-cx", "--center_x", type=float, default=None, required=True,
        help="Coordinate x of box center.")
    parser_dock.add_argument("-cy", "--center_y", type=float, default=None, required=True,
        help="Coordinate y of box center.")
    parser_dock.add_argument("-cz", "--center_z", type=float, default=None, required=True,
        help="Coordinate z of box center.")
    parser_dock.add_argument("-sx", "--size_x", type=float, default=None, required=True,
        help="box size on x.")
    parser_dock.add_argument("-sy", "--size_y", type=float, default=None, required=True,
        help="box size on y.")
    parser_dock.add_argument("-sz", "--size_z", type=float, default=None, required=True,
        help="box size on z.")
    parser_dock.add_argument("--spacing", type=float, default=0.375,
        help="params of ad4")
    parser_dock.add_argument("--ad4_gpf", type=str, default=None,
        help="Grid parameter file for AutoDock4.")
    parser_dock.add_argument("--ad4_dpf", type=str, default=None,
        help="Dock parameter file for AutoDock4.")
    parser_dock.add_argument("-ex", "--exhausiveness", type=int, default=1,
        help="Autodock vina parameter exhausiveness.")
    parser_dock.add_argument("-nm", "--number_of_modes", type=int, default=8,
        help="Autodock vina parameter number_of_modes.")
    parser_dock.add_argument("-er", "--energy_range", type=int, default=3,
        help="Autodock vina parameter energy_range.")
    parser_dock.add_argument("-nga", "--number_of_ga", type=int, default=10,
        help="Autodock4 parameter number_of_ga.")
    parser_dock.add_argument("-maxe", "--maximum_number_of_evals", type=int, default=2500000,
        help="Autodock4 parameter maximum_number_of_evals.")
    parser_dock.add_argument("-maxg", "--maximum_number_of_generations", type=int, default=27000,
        help="Autodock4 parameter maximum_number_of_generations.")
    parser_dock.add_argument("-kn", "--keep_number", type=int, default=None,
        help="Number kept for docking results.")
    parser_dock.add_argument("--not_save_ligs_only", dest="save_ligs_only", action="store_false",
        help="Only save conformations of ligands.")
    parser_dock.add_argument("-noP", "--not_calc_prop", dest="calc_prop", action="store_false",
        help="Not calculate molecular properties.")
    parser_dock.add_argument("-j", "--jobname", type=str, default="dock",
        help="The name of this job.")
    parser_dock.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    parser_dock.add_argument("--projectid", default=None, 
        help="project ID in table")

    parser_MD = subparsers.add_parser("MD", help="molecular dynamics")
    parser_MD.add_argument("--complex_pdb", type=str, required=True,
        help="Complex file in pdb format. If you want to process several complex files at the same time, \
              please seperate them with ','.")
    parser_MD.add_argument("-ts", "--timestep", type=str, default="0.002",
        help="Timestep of simulation. (ps)")
    parser_MD.add_argument("-s", "--steps", type=str, default="500000",
        help="Steps of simulation.")
    parser_MD.add_argument("-pff", "--protein_forcefield", type=str, default="amber99sb",
        help="Forcefield for protein.",
        choices=["amber99sb","amber99sb-ILDN"])
    parser_MD.add_argument("-lff", "--ligand_forcefield", type=str, default="gaff",
        help="Forcefield for ligand.",
        choices=["gaff","gaff2"])
    parser_MD.add_argument("--box_type", type=str, default="triclinic",
        help="Type of box.",
        choices=["cubic", "dodecahedron", "triclinic"])
    parser_MD.add_argument("--boundary", type=str, default="1.0",
        help="Water buffer Size. (nm)")
    parser_MD.add_argument("--heavyh", action="store_true",
        help="Heavy hydrogen atoms (m=4.032)")
    parser_MD.add_argument("--ion_conc", type=str, default="0.15",
        help="Ion concentration. (mol/L)")
    parser_MD.add_argument("-t","--threshold", type=float, default=None,
        help="The RMSD threshold (nm) used to filter the complex.")
    parser_MD.add_argument("-j", "--jobname", type=str, default="MD",
        help="The name of this job.")  
    parser_MD.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    parser_MD.add_argument("-fd", "--from_dock", action="store_true",
        help="Whether use data directly from dock results")
    parser_MD.add_argument("--projectid", default=None, 
        help="project ID in table")

    parser_LBS = subparsers.add_parser("LBS", help="ligand binding site prediction")
    parser_LBS.add_argument("-f", "--protein_file", type=str, required=True,
        help="Protein file in pdb format")
    parser_LBS.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")
    parser_LBS.add_argument("--projectid", default=None, 
        help="project ID in table")

    parser_ADMET = subparsers.add_parser("ADMET", help="ADMET")
    parser_ADMET.add_argument("-f", "--files", required=True, type=str, 
        help="input files")
    parser_ADMET.add_argument("-o", "--output_types", type=str, 
        help="what information we want")
    parser_ADMET.add_argument("--jobname", default="admet", type=str, 
        help="Jobname")
    parser_ADMET.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_ADMET.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_LoopModelling = subparsers.add_parser("LoopModelling", help="loop modelling")
    parser_LoopModelling.add_argument("-p", "--pdbfile", required=True, type=str,
        help="Protein file for loop modelling")
    parser_LoopModelling.add_argument("-s", "--seqfile", default=None, type=str, 
        help="Sequence file for loop modelling. If sequence information is not included in pdb file, you must provide this file!")
    parser_LoopModelling.add_argument("-ft", "--fit_ter", action="store_true",  
        help="whether fit termini")
    parser_LoopModelling.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_LoopModelling.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_LoopOpt = subparsers.add_parser("LoopOpt", help="loop optimization")
    parser_LoopOpt.add_argument("-f", "--pdbfile", required=True, type=str,
        help="Protein file for loop optimization")
    parser_LoopOpt.add_argument("-n_ter", "--start_rnum", required=True, type=int, 
        help="Start residue number.")
    parser_LoopOpt.add_argument("-c_ter", "--end_rnum", required=True, type=int, 
        help="End residue number.")
    parser_LoopOpt.add_argument("-n", "--output_num", default=5, type=int, 
        help="Number of outputs")
    parser_LoopOpt.add_argument("-m", "--method",default="promod", type=str, 
        help="Optimization method")
    parser_LoopOpt.add_argument("-c", "--chain_id", default=None, type=str, 
        help="Chain ID where loop it is if the protein has several chains")
    parser_LoopOpt.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_LoopOpt.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_RID = subparsers.add_parser("RID", help="rid structure refinement")
    parser_RID.add_argument("-f", "--pdbfile", required=True, type=str,
        help="Protein file for refinement")
    parser_RID.add_argument("-nw", "--number_walker", default=4, type=int, 
        help="Number of walkers")
    parser_RID.add_argument("-ni", "--number_iteration", default=2, type=int, 
        help="Number of iterations")
    parser_RID.add_argument("-res", "--selected_residue_index", default="1,2", type=str, 
        help="Residue indice, use ',' to seperate")
    parser_RID.add_argument("-w", "--bottom_width", default=0.6, type=float, 
        help="Bottom width")
    parser_RID.add_argument("-j", "--jobname", default="rid", type=str, 
        help="Jobname")
    parser_RID.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_RID.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_homo = subparsers.add_parser("homo_model", help="homology modelling")
    parser_homo.add_argument("-s", "--seqfile", required=True, type=str,
        help="Protein file for refinement")
    parser_homo.add_argument("-t", "--tmpfile", default=None, type=str,
        help="Protein file for refinement")
    parser_homo.add_argument("-j", "--jobname", default="homology_modelling", type=str, 
        help="Jobname")
    parser_homo.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_homo.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_simsearch = subparsers.add_parser("SimSearch", help="Similarity Search")
    parser_simsearch.add_argument("-s", "--SMILES", required=True, type=str,
        help="smiles string")
    parser_simsearch.add_argument("-d", "--database", default="Targetmol", type=str, 
        choices=["Targetmol", "IBS", "Chembridge", "Enamine", "Bionet", "Specs", "Chemdiv", "Maybridge", "Lifechemicals", "Vitas-m"], 
        help="database name")
    parser_simsearch.add_argument("--scf_m", action="store_false",
        help="scf_m")
    parser_simsearch.add_argument("-j", "--jobname", default="similar", type=str, 
        help="Jobname")
    parser_simsearch.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_simsearch.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_VSW = subparsers.add_parser("VSW", help="vitual screen workflow")
    parser_VSW.add_argument("-p", "--protein_path", required=True, type=str, 
        help="input protein path")
    parser_VSW.add_argument("--jobname", default="vsw", type=str, 
        help="Jobname")
    parser_VSW.add_argument("--projectid", default=None, 
        help="project ID in table")
    parser_VSW.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    parser_query = subparsers.add_parser("QUERY", help="search input param or result oss path")
    parser_query.add_argument("-d", "--type", required=True, type=str, 
        choices=["USER_JOBS", "USER_FILES", "RES_OSS", "PARAM"], help="type of search")
    parser_query.add_argument("--id", default=None, type=int, 
        help="job id for search")
    parser_query.add_argument("-y", "--pass_confirm", action="store_true", 
        help="whether confirm login")

    params = parser.parse_args().__dict__
    if params["cmd"] == "login":
        config_info = {
            "email": params.get("email", None),
            "password": params.get("password", None)}
        sdk_path = os.path.split(os.path.realpath(__file__))[0]
        config_path = os.path.join(sdk_path, "hermite_sdk_config.json")
        with open(config_path, "w") as f:
            json.dump(config_info, f, indent=2)
        hermite = hermite_login(params.get("pass_confirm", False))
        return
    elif params["cmd"] == "project":
        res = createProject(params.get("project", None), params.get("group", None), params.get("pass_confirm", False))
        print(res)
    elif params["cmd"] == "rerun":
        res = rerunTask(params.get("id", None), params.get("pass_confirm", False))
    elif params["cmd"] == "submit":
        submit = submit_file(params)
        submit.upload_file()
    elif params["cmd"] == "proprep":
        protprep = protein_preparation(params)
        protprep.upload_protein()
        protprep.select_group()
        protprep.repair()
        # protprep.protein_prep()
        if params["set_flexible"]:
            protprep.set_flexible()
    elif params["cmd"] == "ligprep":
        ligprep = ligand_preparation(params)
        ligprep.upload_ligand()
        ligprep.ligand_prep()
    elif params["cmd"] == "dock":
        dock = submit_dock(params)
        dock.upload_config()
        dock.submit_dock()
    elif params["cmd"] == "MD":
        md = MD(params)
        md.upload_protein()
        md.submitMD()
    elif params["cmd"] == "LBS":
        lbs = ligand_binding_site_prediction(params)
        lbs.upload_protein()
        lbs.submitLBS()
    elif params["cmd"] == "ADMET":
        admet = ADMET(params)
        admet.upload_files()
        admet.submitADMET()
    elif params["cmd"] == "LoopModelling":
        loop_modelling = LoopModelling(params)
        loop_modelling.upload_file()
        loop_modelling.submitLoopModelling()
    elif params["cmd"] == "LoopOpt":
        loop_opt = LoopOpt(params)
        loop_opt.upload_protein()
        loop_opt.submitLoopOpt()
    elif params["cmd"] == "RID":
        strucure_refinement = StructureRefinement(params)
        strucure_refinement.upload_protein()
        strucure_refinement.submitStructureRefinement()
    elif params["cmd"] == "homo_model":
        homo_modelling = HomologyModelling(params)
        homo_modelling.upload_sequence()
        if params.get("tmpfile", None):
            homo_modelling.upload_template()
        homo_modelling.submitHomologyModelling()
    elif params["cmd"] == "SimSearch":
        simsearch = SimilaritySearch(params)
        simsearch.uploadSimSearch()
    elif params["cmd"] == "VSW":
        vsw = VSW(params)
        vsw.selectTasks()
        vsw.upload_protein(params["protein_path"])
        vsw.uploadParams()
        vsw.submitVSW()
    elif params["cmd"] == "QUERY":
        query = QUERY(params)
        query.query()
    return


if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    print("Time: ", "{}s".format(end_time - start_time))
