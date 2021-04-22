import sys
from BTrees.OOBTree import OOBTree
from suffix_trees import STree

class stdict(OOBTree):
    """
    Custom class to allow substring matches in keys.
    """
    def __init__(self, *args, **kwargs):
        super(stdict, self).__init__(*args, **kwargs)
        self.key = ''
        self.stree = STree.STree([])
        self.indices = {}

    def __setitem__(self, key, value):
        if key:
            inKeys = self.stree.find_all(key)
            if inKeys and self.indices:
                for k in inKeys:

                    start = super(stdict, self).maxKey(k)
                    end = self.indices[start]

                    if self.key[start:self.indices[start] + 1] == key:
                        raise KeyError("Error, key already in stdict")

            start = len(self.key)
            if not start == 0:
                end = start + len(key) - 1
            else:
                end = 0
            self.key += key
            self.stree = STree.STree(self.key)
            self.indices[start] = end

            super(stdict, self).__setitem__(start, value)

    def __getitem__(self, approxkey):
        if not approxkey:
            raise KeyError(approxkey)

        idxs = self.stree.find_all(approxkey)

        if not idxs:
            raise KeyError(approxkey)

        else:
            startIdxs = set(super(stdict, self).maxKey(i) for i in idxs)
            ret = []
            for sid in startIdxs:
                item = super(stdict, self).__getitem__(sid)
                if not item is None:
                    ret.append(item)

            return ret

            #return [super(stdict, self).__getitem__(sid) for sid in startIdxs]

    def __call__(self, key, value):
        idxs = self.stree.find_all(key)
        if not idxs:
            raise KeyError(key)

        pair = tuple()
        minLen = sys.maxint
        for i in idxs:
            startIdx = super(stdict, self).maxKey(i)
            endIdx = self.indices[startIdx]
            pLen = endIdx - startIdx

            if pLen < minLen:
                pair = (startIdx, endIdx)
                minLen = pLen

        super(stdict, self).__setitem__(pair[0], None)
        replace = (pair[1] - pair[0] + 1)*'X'
        self.key =  ''.join((self.key[:pair[0]], replace, self.key[pair[1] + 1:]))
        start = len(self.key)
        end = start + len(key) - 1
        self.key += key
        self.stree = STree.STree(self.key)
        self.indices[start] = end

        super(stdict, self).__setitem__(start, value)

    def __iter__(self):
        sk = memoryview(self.key)
        for key in super(stdict, self).keys():
            rKey = sk[key:self.indices[key] + 1].tobytes()
            if 'X' not in rKey:
                yield sk[key:self.indices[key] + 1].tobytes()