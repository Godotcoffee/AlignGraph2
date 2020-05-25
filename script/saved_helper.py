import os

def check(root_dir, *file_list):
    check_path = os.path.join(root_dir, 'CHECK')
    if not os.path.isfile(check_path):
        return False
    file_dict = dict()
    try:
        with open(check_path) as check:
            for line in check:
                sp = line.strip('\n').split('\t')
                file_dict[sp[0]] = (sp[1], sp[2])
        for check_file in file_list:
            abs_path = os.path.abspath(check_file)
            if abs_path not in file_dict:
                return False
            t = file_dict[abs_path]
            if t[0] != str(os.path.getsize(abs_path)) or t[1] != str(os.path.getmtime(abs_path)):
                return False
    except:
        return False
    return True

def save(root_dir, *file_list):
    check_path = os.path.join(root_dir, 'CHECK')
    try:
        with open(check_path, 'w') as check:
            for check_file in file_list:
                abs_path = os.path.abspath(check_file)
                check.write('{}\t{}\t{}\n'.format(abs_path, str(os.path.getsize(abs_path)), str(os.path.getmtime(abs_path))))
    except:
        pass


def remove(root_dir):
    check_path = os.path.join(root_dir, 'CHECK')
    try:
        if os.path.exists(check_path):
            os.remove(check_path)
    except:
        pass


def check_args(root_dir, **kw):
    check_args_path = os.path.join(root_dir, 'ARGS')
    if not os.path.isfile(check_args_path):
        return False
    args_dict = dict()
    try:
        with open(check_args_path) as check:
            for line in check:
                sp = line.strip('\n').split('\t')
                args_dict[sp[0]] = sp[1]
        for k in kw:
            if k not in args_dict:
                return False
            if args_dict[k] != str(kw[k]):
                return False
    except:
        return False
    return True

def save_args(root_dir, **kw):
    check_args_path = os.path.join(root_dir, 'ARGS')
    try:
        with open(check_args_path, 'w') as check:
            for k in kw:
                check.write('{}\t{}\n'.format(k, kw[k]))
    except:
        pass


def remove_args(root_dir):
    check_args_path = os.path.join(root_dir, 'ARGS')
    try:
        if os.path.exists(check_args_path):
            os.remove(check_args_path)
    except:
        pass


if __name__ == "__main__":
    print(check_args('.', k=1, kk='23', o=3.2))
    save_args('.', k=1, kk='23', o=3.2)
