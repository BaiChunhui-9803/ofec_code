#include"hash_generator.h"
#include"crc32.h"
#include"keccak.h"
#include"md5.h"
#include"sha1.h"
#include"sha256.h"
#include"sha3.h"

void HashGenerator(std::unique_ptr<Hash>& hash_fun, HashType type) {
	switch (type)
	{
	case HashType::CRC32: {
		hash_fun.reset(new CRC32());
		break;
	}
	case HashType::Keccak: {
		hash_fun.reset(new Keccak());
		break;
	}
	case HashType::MD5: {
		hash_fun.reset(new MD5());
		break;
	}
	case HashType::SHA1: {
		hash_fun.reset(new SHA1());
		break;
	}
	case HashType::SHA256: {
		hash_fun.reset(new SHA256());
		break;
	}
	case HashType::SHA3: {
		hash_fun.reset(new SHA3());
		break;
	}
	default:
		hash_fun.reset(nullptr);
		break;
	}
}


