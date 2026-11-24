#!/usr/bin/env python3
"""Streaming reader for cil_bd / cil_p_bd snapshots (little endian).
CILTRJ03 (current): float positions relative to each frame's centre of mass (returned as
unwrapped doubles, ~1e-6 absolute), bonds as full sorted lists on keyframes and as changes
otherwise (returned as the full sorted list every frame), and for CIL-P the rest lengths of
bonds not at sigma_ij (returned as frame['rest_lengths'] for every bond).
Older formats: CILTRJ01 (CIL-C) and CILTRJ02 (CIL-P, rest length of every bond).
Flags: 1 structural frame, 2 burst frame (a frame may be both), 4 keyframe.
Burst origins and frame step are relative to the start of measurement.
Pair indices are mapped to particle IDs by bonds_to_pairs(). Read frames in order.
"""
import argparse,struct
from pathlib import Path
import numpy as np

def exact(f,n):
    b=f.read(n)
    if len(b)!=n:raise ValueError('Truncated trajectory')
    return b

def bonds_to_pairs(indices,n):
    indices=np.asarray(indices,dtype=np.int64)
    if np.any(indices<0) or np.any(indices>=n*(n-1)//2):raise ValueError('Invalid bond index')
    i=np.arange(n,dtype=np.int64);starts=i*n-i*(i+1)//2
    a=np.searchsorted(starts,indices,side='right')-1
    return np.column_stack([a,a+1+indices-starts[a]])

def read_frames(path):
    with open(path,'rb') as f:
        magic=exact(f,8)
        if magic not in (b'CILTRJ01',b'CILTRJ02',b'CILTRJ03'):raise ValueError('Unknown trajectory format')
        endian,n=struct.unpack('<II',exact(f,8))
        if endian!=0x01020304 or not 0<n<=16384:raise ValueError('Invalid trajectory header')
        model=struct.unpack('<I',exact(f,4))[0] if magic==b'CILTRJ03' else int(magic==b'CILTRJ02')
        radii=np.frombuffer(exact(f,8*n),dtype='<f8').copy()
        previous=-1;bonds=None
        while True:
            if magic==b'CILTRJ03':
                raw=f.read(56)
                if not raw:break
                if len(raw)!=56:raise ValueError('Partial frame header')
                step,t,lx,ly,flags,count,nb,nrem,nadd,nr0=struct.unpack('<qdddIIIIII',raw)
            else:
                raw=f.read(44)
                if not raw:break
                if len(raw)!=44:raise ValueError('Partial frame header')
                step,t,lx,ly,flags,count,nb=struct.unpack('<qdddIII',raw)
            if step<=previous or count>16 or nb>n*(n-1)//2:raise ValueError('Invalid frame header')
            previous=step
            origins=np.frombuffer(exact(f,8*count),dtype='<i8').copy()
            if magic==b'CILTRJ03':
                cx,cy=struct.unpack('<dd',exact(f,16))
                xy=np.frombuffer(exact(f,8*n),dtype='<f4').reshape(n,2).astype(np.float64)+(cx,cy)
                if flags&4:bonds=np.frombuffer(exact(f,4*nb),dtype='<u4').copy()
                else:
                    if bonds is None:raise ValueError('Bond changes before the first keyframe')
                    rem=np.frombuffer(exact(f,4*nrem),dtype='<u4');add=np.frombuffer(exact(f,4*nadd),dtype='<u4')
                    bonds=np.union1d(np.setdiff1d(bonds,rem,assume_unique=True),add).astype(np.uint32)
                if len(bonds)!=nb:raise ValueError('Bond changes do not add up')
                over=np.frombuffer(exact(f,8*nr0),dtype=[('pair','<u4'),('r0','<f4')])
            else:
                xy=np.frombuffer(exact(f,16*n),dtype='<f8').reshape(n,2).copy()
                bonds=np.frombuffer(exact(f,4*nb),dtype='<u4').copy()
            frame=dict(step=step,time=t,box=(lx,ly),flags=flags,origins=origins,
                       positions=xy,radii=radii,bonds=bonds)
            if magic==b'CILTRJ02':frame['rest_lengths']=np.frombuffer(exact(f,8*nb),dtype='<f8').copy()
            elif model==1:
                p=bonds_to_pairs(bonds,n);r0=radii[p[:,0]]+radii[p[:,1]]
                if len(over):r0[np.searchsorted(bonds,over['pair'])]=over['r0']
                frame['rest_lengths']=r0
            yield frame

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('file',type=Path);args=ap.parse_args()
    count=sparse=burst=0;first=last=None
    for frame in read_frames(args.file):
        count+=1;sparse+=bool(frame['flags']&1);burst+=bool(frame['flags']&2)
        if first is None:first=frame['time']
        last=frame['time'];bonds_to_pairs(frame['bonds'],len(frame['radii']))
    print(f'{count} frames; {sparse} structural, {burst} burst; time {first} .. {last}; {args.file.stat().st_size/1e6:.3f} MB')
