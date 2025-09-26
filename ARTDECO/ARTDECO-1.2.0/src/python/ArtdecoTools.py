
def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

def is_int(value):
    try:
        if int(float(value)) == float(value):
            return True
        else:
            return False
    except ValueError:
        return False

def clean_list(liste):
    result = []
    for i in range(len(liste)):
        if is_float(liste[i]):
            result.append(liste[i])
    return result

def clean_line(line):
    result = []
    liste = line.split(None)
    for item in liste:
        if is_float(item):
            result.append(item)
    return result

def clean_comment_lines(liste):
    result = []
    for line in liste:
        if line[0] != '#':
            result.append(line)
    return result
