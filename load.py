import gmpy2
import re
from structure import AffinePointG1,AffinePointG2,G2Coordinate,UniversalParams
from bls12_381 import fq,fr
from plonk_core.src.proof_system.prover_key import Prover_Key
from plonk_core.src.proof_system.widget.arithmetic import Arith
from plonk_core.src.proof_system.widget.lookup import Lookup
from plonk_core.src.proof_system.permutation import Permutation
def parse_bigint(s):
    start = s.find('"(') + 2
    end = s.find(')')
    bigint_str = s[start:end]
    return gmpy2.mpz(bigint_str,16)

def read_pp_data(filename):
    # 打开文本文件以读取数据
    with open(filename, "r") as file:
        data = file.read()

    powers_of_g = []
    powers_of_gamma_g = []
    h=AffinePointG2(x=G2Coordinate(c0=gmpy2.mpz(1), c1=gmpy2.mpz(2)), y=G2Coordinate(c0=gmpy2.mpz(3), c1=gmpy2.mpz(4)))
    beta_h=AffinePointG2(x=G2Coordinate(c0=gmpy2.mpz(1), c1=gmpy2.mpz(2)), y=G2Coordinate(c0=gmpy2.mpz(3), c1=gmpy2.mpz(4)))

    lines = data.split('\n')
    current_section = None

    for line in lines:
        if line.endswith(":"):
            current_section = line.rstrip(":")
        elif line.startswith("["):
            values = line.strip("[] ").split()
            if current_section == "powers_of_g":
                x_str = values[1]
                y_str = values[2]
                G1_point= AffinePointG1(x=fq.Fq(value=parse_bigint(x_str)),y=fq.Fq(value=parse_bigint(y_str)))
                powers_of_g.append(G1_point)
            elif current_section == "powers_of_gamma_g":
                x_str = values[1]
                y_str = values[2]
                G1_point= AffinePointG1(x=fq.Fq(value=parse_bigint(x_str)),y=fq.Fq(value=parse_bigint(y_str)))
                powers_of_gamma_g.append(G1_point)
            elif current_section == "h":
                h.x.c0 = fq.Fq(value=parse_bigint(values[1]))
                h.x.c1 = fq.Fq(value=parse_bigint(values[2]))
                h.y.c0 = fq.Fq(value=parse_bigint(values[3]))
                h.y.c1 = fq.Fq(value=parse_bigint(values[4]))
                
            elif current_section == "beta_h":
                beta_h.x.c0 = fq.Fq(value=parse_bigint(values[1]))
                beta_h.x.c1 = fq.Fq(value=parse_bigint(values[2]))
                beta_h.y.c0 = fq.Fq(value=parse_bigint(values[3]))
                beta_h.y.c1 = fq.Fq(value=parse_bigint(values[4]))
                
    pp = UniversalParams(powers_of_g, powers_of_gamma_g, h, beta_h)
    return pp

def read_pk_data(filename):
    
    with open(filename, "r") as file:
        lines = file.readlines()

    n = int(lines[0].split(':')[1].strip())
    data = {}
    data["n"]=n
    current_key = None

    # 解析文件数据
    for line in lines[1:]:
        line = line.strip()
        if line.startswith("arithmetic:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("range_selector:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("logic_selector:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("lookup:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("fixed_group_add_selector:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("variable_group_add_selector:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("permutation:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("linear_evaluations:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("v_h_coset_8n:"):
            current_key = line[:-1]
            data[current_key] = {}
            subkey=None
        elif line.startswith("q_m"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_l"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_r"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_o"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_4"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_c"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_hl"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_hr"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_h4"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("q_arith"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.endswith(":") and current_key=="range_selector" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="range_selector" and subkey:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="logic_selector" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="logic_selector" and subkey:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.startswith("q_lookup:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("table1:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("table2:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("table3:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("table4:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("left_sigma:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("right_sigma:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("out_sigma:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.startswith("fourth_sigma:"):
            subkey=line[:-1]
            data[current_key][subkey]={}
            subsubkey=None
        elif line.endswith(":") and current_key=="fixed_group_add_selector" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="fixed_group_add_selector" and subkey:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="variable_group_add_selector" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="variable_group_add_selector" and subkey:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="linear_evaluations" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key=="v_h_coset_8n" and subkey==None:
            subkey=line[:-1]
            data[current_key][subkey]=[]
            subsubkey=None
        elif line.endswith(":") and current_key and subkey and subsubkey==None:
            subsubkey=line[:-1]
            data[current_key][subkey][subsubkey]=[]
        elif line.endswith(":") and current_key and subkey and subsubkey:
            subsubkey=line[:-1]
            data[current_key][subkey][subsubkey]=[]
        elif not(line.endswith(":")):
            value = parse_bigint(line)
            value = fr.Fr.from_repr(value)
            if subkey==None and subsubkey==None:
                data[current_key].append(value)
            elif subkey and subsubkey==None:
                data[current_key][subkey].append(value)
            elif line.startswith("-") and subsubkey:
                subsubkey="coeffs"
                data[current_key][subkey][subsubkey].append(value)
            elif subsubkey:
                subsubkey="evals"
                data[current_key][subkey][subsubkey].append(value)
    
    arithmetic = Arith(q_m=(data["arithmetic"]["q_m"]["coeffs"],data["arithmetic"]["q_m"]["evals"]),
                       q_l=(data["arithmetic"]["q_l"]["coeffs"],data["arithmetic"]["q_l"]["evals"]),
                       q_r=(data["arithmetic"]["q_r"]["coeffs"],data["arithmetic"]["q_r"]["evals"]),
                       q_o=(data["arithmetic"]["q_o"]["coeffs"],data["arithmetic"]["q_o"]["evals"]),
                       q_4=(data["arithmetic"]["q_4"]["coeffs"],data["arithmetic"]["q_4"]["evals"]),
                       q_c=(data["arithmetic"]["q_c"]["coeffs"],data["arithmetic"]["q_c"]["evals"]),
                       q_hl=(data["arithmetic"]["q_hl"]["coeffs"],data["arithmetic"]["q_hl"]["evals"]),
                       q_hr=(data["arithmetic"]["q_hr"]["coeffs"],data["arithmetic"]["q_hr"]["evals"]),
                       q_h4=(data["arithmetic"]["q_h4"]["coeffs"],data["arithmetic"]["q_h4"]["evals"]),
                       q_arith=(data["arithmetic"]["q_arith"]["coeffs"],data["arithmetic"]["q_arith"]["evals"]))
    
    lookup = Lookup(q_lookup=(data["lookup"]["q_lookup"]["coeffs"],data["lookup"]["q_lookup"]["evals"]),
                    table_1=data["lookup"]["table1"]["coeffs"],
                    table_2=data["lookup"]["table2"]["coeffs"],
                    table_3=data["lookup"]["table3"]["coeffs"],
                    table_4=data["lookup"]["table4"]["coeffs"])
    
    permutation = Permutation(left_sigma=(data["permutation"]["left_sigma"]["coeffs"],data["permutation"]["left_sigma"]["evals"]),
                              right_sigma=(data["permutation"]["right_sigma"]["coeffs"],data["permutation"]["right_sigma"]["evals"]),
                              out_sigma=(data["permutation"]["out_sigma"]["coeffs"],data["permutation"]["out_sigma"]["evals"]),
                              fourth_sigma=(data["permutation"]["fourth_sigma"]["coeffs"],data["permutation"]["fourth_sigma"]["evals"]),
                              linear_evaluations=data["linear_evaluations"]["evals"])
    
    pk = Prover_Key(arithmetic=arithmetic,range_selector=(data["range_selector"]["coeffs"],data["range_selector"]["evals"]),
                    logic_selector=(data["logic_selector"]["coeffs"],data["logic_selector"]["evals"]),lookup=lookup,
                    fixed_group_add_selector=(data["fixed_group_add_selector"]["coeffs"],data["fixed_group_add_selector"]["evals"]),
                    variable_group_add_selector=(data["variable_group_add_selector"]["coeffs"],data["variable_group_add_selector"]["evals"]),
                    permutation=permutation,v_h_coset_8n=data["v_h_coset_8n"]["evals"])
    
    return pk

def read_cs_data(filename):
    with open(filename,"r")as file:
        lines = file.readlines()
    n = int(lines[0].split(':')[1].strip())
    intended_pi_pos=lines[1].split(':')[1].strip()
    zero_var=int(lines[2].split(':')[1].strip())
    pattern = r'\[([\d\s,]+)\]'
    match = re.search(pattern, lines[-1])
    pi = match.group(1)
    public_inputs_list = [gmpy2.mpz(x) for x in pi.split(',')]
    public_inputs = (public_inputs_list[3] << 192) | (public_inputs_list[2] << 128) | (public_inputs_list[1] << 64) | public_inputs_list[0]
    public_inputs = fr.Fr(value=public_inputs)
    data = {}
    data["n"]=n
    data["intended_pi_pos"]=eval(intended_pi_pos)
    data["zero_var"]=zero_var
    data["public_inputs"]=public_inputs
    current_key = None

    # 解析文件数据
    for line in lines[3:-2]:
        line = line.strip()
        if line.endswith(":"):
            current_key=line.rstrip(":")
            data[current_key]=[]
        elif line.startswith("Fp"):
            value = parse_bigint(line)
            value = fr.Fr(value=value)
            data[current_key].append(value)
        else:
            data[current_key].append(int(line))
    return data

def read_scalar_data(filename):
    with open(filename, "r") as file:
        content = file.read()

    # 使用正则表达式匹配包含列表的字符串
    pattern = r'\[.*?\]'
    matches = re.findall(pattern, content)

    big_list = []

    # 如果第一个元素不符合正则表达式，手动添加到大列表
    if content.startswith('Fp256'):
        str_list=content.split('Fp256(BigInteger256(')[1].split(')')[0]
        element_str = str_list.strip('[]')
        elements = [gmpy2.mpz(e) for e in element_str.split(',')]
        combined_value = 0
        for u64_element in reversed(elements):
            combined_value = (combined_value << 64) | u64_element
        big_list.append(fr.Fr(value=gmpy2.mpz(combined_value)))

    matches.pop(0)

    # 遍历匹配的字符串列表并处理
    for match in matches:
        # 去除字符串中的方括号，并按逗号分割成元素列表
        elements_str = match.strip('[]')
        elements = [gmpy2.mpz(e) for e in elements_str.split(',')]
        # 将元素列表作为子列表追加到大列表中
        # 将元素列表合并成一个256位整数
        combined_value = 0
        for u64_element in reversed(elements):
            combined_value = (combined_value << 64) | u64_element
        big_list.append(fr.Fr(value=gmpy2.mpz(combined_value)))
    return big_list