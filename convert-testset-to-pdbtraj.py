#!/usr/bin/env python
#
# Quick and dirty script to convert Itay Barel's lipid membrane test set into
# something that VMD can load, in order to test bending modulus plugin
#
# Tom Joseph <thomas.joseph@uphs.upenn.edu>
import argparse

def main():
	ap = argparse.ArgumentParser(description='Just read the source code')
	ap.add_argument('num_lipids', type=int)
	# ap.add_argument('num_frames', type=int)
	ap.add_argument('head_name', nargs='?', type=str, default='P')
	ap.add_argument('tail_name', nargs='?', type=str, default='C218')
	args = ap.parse_args()

	boxx = open('boxsizeX.out').readlines()
	boxy = open('boxsizeY.out').readlines()
	boxz = open('boxsizeZ.out').readlines()
	lipidx = open('LipidX.out').readlines()
	lipidy = open('LipidY.out').readlines()
	lipidz = open('LipidZ.out').readlines()

	print('Lines in box files: {} {} {}'.format(len(boxx), len(boxy), len(boxz)))

	if len(boxx) != len(boxy) or len(boxy) != len(boxz):
		print('Lengths of box files are not equal. I am confused.')
		exit(1)

	if len(lipidx) != len(lipidy) or len(lipidy) != len(lipidz):
		print('Lengths of lipid files are not equal. I am confused.')
		exit(1)

	print('Lines in lipid files: {} {} {}'.format(len(lipidx), len(lipidy), len(lipidz)))
	num_frames = len(lipidx) / args.num_lipids / 2.0
	if int(num_frames) != num_frames:
		print('Non-integral number of frames: {}'.format(num_frames))
		exit(1)

	num_frames = int(num_frames)
	print('Number of frames is therefore: {}'.format(num_frames))
	if len(boxx) != num_frames:
		print('But the number of box lines is not the same as the calculated number of frames. Tread carefully.')

	print('Lipid heads will be named {} and tails will be named {}.'.format(args.head_name, args.tail_name))

	out = open('lipid-traj.pdb', 'w')
	cryst1 = 'CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}{:10s}{:3d}'
	atomrec = '{:6s}{:5d} {:^4s}{:1s}{:4s}{:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'

	for frame in range(num_frames):
		# CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1
		a, b, c = float(boxx[frame]), float(boxy[frame]), float(boxz[frame])
		out.write(cryst1.format(a, b, c, 90.0, 90.0, 90.0, 'P', 1))
		out.write('\n')


		for lipid in range(args.num_lipids):
			offset = frame * args.num_lipids * 2 + lipid * 2
			hx, hy, hz = float(lipidx[offset]), float(lipidy[offset]), float(lipidz[offset])
			tx, ty, tz = float(lipidx[offset+1]), float(lipidy[offset+1]), float(lipidz[offset+1])

			atomid = lipid * 2 + 1
			resid = lipid + 1
			lipid_name = 'POPC'

            # "ATOM", self.curr_atom_id, atom.name, atom.altLoc, atom.residue.resname, atom.residue.segid, resid,
            #       ' ', atom.position[0], atom.position[1], atom.position[2],
            #       occupancy, tempfactor, atom.type, '  ')) # occupancy, temp factor, type, charge
			out.write(atomrec.format('ATOM', atomid, args.head_name, ' ', lipid_name, 'X', resid, ' ', hx, hy, hz,
				1.0, 0.0, args.head_name[0], ' '))

			out.write(atomrec.format('ATOM', atomid, args.tail_name, ' ', lipid_name, 'X', resid, ' ', tx, ty, tz,
				1.0, 0.0, args.tail_name[0], ' '))

		out.write('TER\nEND\n')

if __name__ == '__main__':
	main()