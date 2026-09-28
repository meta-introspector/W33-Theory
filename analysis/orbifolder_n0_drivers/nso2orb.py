"""convert a non-SUSY-orbifolder (Z2W x Z_M) geometry/model to orbifolder-1.2.1 Z_M x Z_N format (drop the unused Z_K index)"""
import sys, re
def geom(src, dst):
    L = open(src).read().splitlines()
    out = []; i = 0
    while i < len(L):
        l = L[i]
        if l.strip() == 'point group':
            out += [l, L[i+1], L[i+2]]; i += 4; continue           # drop K
        if l.strip() == 'ZMxZNxZK':
            i += 2; continue
        f = l.split()
        if len(f) == 9 and all(re.fullmatch(r'-?\d+', x) for x in f):
            assert f[2] == '0', l
            out.append(' '.join(f[:2] + f[3:])); i += 1; continue
        out.append(l); i += 1
    open(dst, 'w').write('\n'.join(out) + '\n')
def model(src, dst, geomname):
    t = open(src).read()
    blocks = []
    for blk in t.split('begin model')[1:]:
        head, rest = blk.split('Shifts and Wilsonlines:')
        rows = rest.split('end model')[0].strip().splitlines()
        assert len(rows) == 9 and all(x.strip(',') in ('0/1', '0') for x in rows[2].split()), rows[2]
        rows = rows[:2] + rows[3:]
        head = re.sub(r'SpaceGroup:\S+', 'SpaceGroup:Geometry/' + geomname, head)
        blocks.append('begin model' + head + 'Shifts and Wilsonlines:\n' + '\n'.join(rows) + '\nend model\n')
    open(dst, 'w').write('\n'.join(blocks))
if __name__ == '__main__':
    if sys.argv[1] == 'geom': geom(sys.argv[2], sys.argv[3])
    else: model(sys.argv[2], sys.argv[3], sys.argv[4])
