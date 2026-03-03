#!/usr/bin/env python3
import subprocess

def run_tmalign(query_pdb, template_pdb):
    """
    TM-align을 실행하여 stdout 문자열을 반환
    """
    cmd = ["/applic/bin/TMalign", query_pdb, template_pdb]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout


def parse_tmalign_matrix(text):
    """
    TM-align stdout에서 회전행렬 + translation 추출
    """
    lines = text.splitlines()
    t = [0.0, 0.0, 0.0]
    u = [[0.0]*3 for _ in range(3)]

    idx = None
    for i, line in enumerate(lines):
        if "Rotation matrix to rotate Chain-1 to Chain-2" in line:
            idx = i
            break
    if idx is None:
        raise RuntimeError("Rotation matrix not found")

    # 실제 matrix는 idx+2 ~ idx+4
    for i in range(3):
        parts = lines[idx+2+i].split()
        t[i] = float(parts[1])
        u[i][0] = float(parts[2])
        u[i][1] = float(parts[3])
        u[i][2] = float(parts[4])

    return t, u


def transform_coord(x, y, z, t, u):
    X = t[0] + u[0][0]*x + u[0][1]*y + u[0][2]*z
    Y = t[1] + u[1][0]*x + u[1][1]*y + u[1][2]*z
    Z = t[2] + u[2][0]*x + u[2][1]*y + u[2][2]*z
    return X, Y, Z


def apply_transform_to_pdb(in_pdb, out_pdb, t, u):
    with open(in_pdb, "r") as fin, open(out_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except:
                    fout.write(line)
                    continue
                X, Y, Z = transform_coord(x, y, z, t, u)

                new_line = (
                    line[:30]
                    + f"{X:8.3f}{Y:8.3f}{Z:8.3f}"
                    + line[54:]
                )
                fout.write(new_line)
            else:
                fout.write(line)


def align_pdb(query, template, out_pdb):
    # 1) TM-align 실행
    stdout = run_tmalign(query, template)

    # 2) matrix 추출
    t, u = parse_tmalign_matrix(stdout)

    # 3) query 좌표를 template 좌표계로 변환
    apply_transform_to_pdb(query, out_pdb, t, u)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python align.py <query.pdb> <template.pdb> <output.pdb>")
        sys.exit(1)

    query = sys.argv[1]
    template = sys.argv[2]
    out_pdb = sys.argv[3]

    align_pdb(query, template, out_pdb)
