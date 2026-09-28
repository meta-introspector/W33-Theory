"""convert an orbifolder-1.2.1 Z_N geometry/model into the non-SUSY-orbifolder format (point group M,N,K; 9-column rows)"""
import sys, re
def geom(src, dst):
    L = open(src).read().splitlines()
    out = []; i = 0
    while i < len(L):
        l = L[i]
        if l.strip() == 'begin discrete symmetries':
            while L[i].strip() != 'end discrete symmetries': i += 1
            i += 1; continue
        if l.strip() == 'point group':
            out += [l, L[i+1], L[i+2] if L[i+2].strip() else '1', '1']
            i += 3 if L[i+2].strip() else 2
            continue
        if l.strip() == 'ZMxZN':
            out += [l, L[i+1], '', 'ZMxZNxZK', 'false']; i += 2; continue
        f = l.split()
        if len(f) == 8 and all(re.fullmatch(r'-?\d+', x) for x in f):
            out.append(' '.join(f[:2] + ['0'] + f[2:])); i += 1; continue
        out.append(l); i += 1
    open(dst, 'w').write('\n'.join(out) + '\n')
def model(src, dst, geomname):
    t = open(src).read()
    blocks = []
    for blk in t.split('begin model')[1:]:
        head, rest = blk.split('Shifts and Wilsonlines:')
        rows = rest.split('end model')[0].strip().splitlines()
        assert len(rows) == 8
        rows = rows[:2] + [' '.join(['0'] * 16)] + rows[2:]
        head = re.sub(r'SpaceGroup:\S+', 'SpaceGroup:Geometry/' + geomname, head)
        blocks.append('begin model' + head + 'Shifts and Wilsonlines:\n' + '\n'.join(rows) + '\nend model\n')
    open(dst, 'w').write('\n'.join(blocks))
if __name__ == '__main__':
    if sys.argv[1] == 'geom': geom(sys.argv[2], sys.argv[3])
    else: model(sys.argv[2], sys.argv[3], sys.argv[4])
