import math
from itertools import product
from torch import Tensor
import torch
import torch.nn as nn
from torch.nn import Transformer
from torchtext.vocab import build_vocab_from_iterator
import torch.nn.functional as F

DEVICE = torch.device("cpu")
EOS_IDX = 3 


# def create_vocab_and_transforms():

#     protein_vocab = list("ACDEFGHIKLMNPQRSTVWY*")
#     codon_vocab = list("".join(x) for x in product("ATCG", repeat=3))

#     SRC_LANGUAGE = "protein"
#     TGT_LANGUAGE = "codons"

#     # Place-holders
#     token_transform = {}
#     vocab_transform = {}

#     UNK_IDX, PAD_IDX, BOS_IDX, EOS_IDX = 0, 1, 2, 3
#     special_symbols = ['<unk>', '<pad>', '<bos>', '<eos>']

#     # first the protein vocab!
#     vocab_transform[SRC_LANGUAGE] = build_vocab_from_iterator([protein_vocab], min_freq=1, specials=special_symbols, special_first=True)

#     # now the codon vocab!
#     vocab_transform[TGT_LANGUAGE] = build_vocab_from_iterator([codon_vocab], min_freq=1, specials=special_symbols, special_first=True)

#     for lang in [SRC_LANGUAGE, TGT_LANGUAGE]:
#         vocab_transform[lang].set_default_index(UNK_IDX)

#     # now the token transforms
#     token_transform[SRC_LANGUAGE] = lambda x: x.split(" ")
#     token_transform[TGT_LANGUAGE] = lambda x: x.split(" ")

#     return vocab_transform, token_transform



# functions for generating the autoregressive mask
def generate_square_subsequent_mask(sz):
    mask = (torch.triu(torch.ones((sz, sz), device=DEVICE)) == 1).transpose(0, 1)
    mask = mask.float().masked_fill(mask == 0, float('-inf')).masked_fill(mask == 1, float(0.0))
    return mask


def create_mask(src, tgt):
    src_seq_len = src.shape[0]
    tgt_seq_len = tgt.shape[0]

    tgt_mask = generate_square_subsequent_mask(tgt_seq_len)
    src_mask = torch.zeros((src_seq_len, src_seq_len),device=DEVICE).type(torch.bool)

    src_padding_mask = (src == PAD_IDX).transpose(0, 1)
    tgt_padding_mask = (tgt == PAD_IDX).transpose(0, 1)
    return src_mask, tgt_mask, src_padding_mask, tgt_padding_mask


# helper Module that adds positional encoding to the token embedding to introduce a notion of word order.
class PositionalEncoding(nn.Module):
    def __init__(self,
                 emb_size: int,
                 dropout: float,
                 maxlen: int = 5000):
        super(PositionalEncoding, self).__init__()
        den = torch.exp(- torch.arange(0, emb_size, 2)* math.log(10000) / emb_size)
        pos = torch.arange(0, maxlen).reshape(maxlen, 1)
        pos_embedding = torch.zeros((maxlen, emb_size))
        pos_embedding[:, 0::2] = torch.sin(pos * den)
        pos_embedding[:, 1::2] = torch.cos(pos * den)
        pos_embedding = pos_embedding.unsqueeze(-2)

        self.dropout = nn.Dropout(dropout)
        self.register_buffer('pos_embedding', pos_embedding)

    def forward(self, token_embedding: Tensor):
        return self.dropout(token_embedding + self.pos_embedding[:token_embedding.size(0), :])


# helper Module to convert tensor of input indices into corresponding tensor of token embeddings
class TokenEmbedding(nn.Module):
    def __init__(self, vocab_size, emb_size):
        super(TokenEmbedding, self).__init__()
        self.embedding = nn.Embedding(vocab_size, emb_size)
        self.emb_size = emb_size

    def forward(self, tokens):
        return self.embedding(tokens.long()) * math.sqrt(self.emb_size)


# Seq2Seq Network
class Seq2SeqTransformer(nn.Module):
    def __init__(self,
                 num_encoder_layers: int,
                 num_decoder_layers: int,
                 emb_size: int,
                 nhead: int,
                 src_vocab_size: int,
                 tgt_vocab_size: int,
                 dim_feedforward: int = 512,
                 dropout: float = 0.1):
        super(Seq2SeqTransformer, self).__init__()
        self.transformer = Transformer(d_model=emb_size,
                                       nhead=nhead,
                                       num_encoder_layers=num_encoder_layers,
                                       num_decoder_layers=num_decoder_layers,
                                       dim_feedforward=dim_feedforward,
                                       dropout=dropout)
        self.generator = nn.Linear(emb_size, tgt_vocab_size)
        self.src_tok_emb = TokenEmbedding(src_vocab_size, emb_size)
        self.tgt_tok_emb = TokenEmbedding(tgt_vocab_size, emb_size)
        self.positional_encoding = PositionalEncoding(
            emb_size, dropout=dropout)

    def forward(self, src, tgt, src_mask, tgt_mask, src_padding_mask, tgt_padding_mask, memory_key_padding_mask):
        src_emb = self.positional_encoding(self.src_tok_emb(src))
        tgt_emb = self.positional_encoding(self.tgt_tok_emb(tgt))
        outs = self.transformer(src_emb, tgt_emb, src_mask, tgt_mask, None, src_padding_mask, tgt_padding_mask, memory_key_padding_mask)
        return self.generator(outs)

    def encode(self, src, src_mask):
        return self.transformer.encoder(self.positional_encoding(self.src_tok_emb(src)), src_mask)

    def decode(self, tgt, memory, tgt_mask):
        return self.transformer.decoder(self.positional_encoding(self.tgt_tok_emb(tgt)), memory, tgt_mask)

    def sample(self, src, src_mask, max_len, start_symbol, temperature=1.0):
        src = src.to(DEVICE)
        src_mask = src_mask.to(DEVICE)

        memory = self.encode(src, src_mask)
        ys = torch.ones(1, 1).fill_(start_symbol).type(torch.long).to(DEVICE)
        for i in range(max_len-1):
            memory = memory.to(DEVICE)
            tgt_mask = (generate_square_subsequent_mask(ys.size(0))
                        .type(torch.bool)).to(DEVICE)
            out = self.decode(ys, memory, tgt_mask)
            out = out.transpose(0, 1)
            logits = self.generator(out[:, -1]) / temperature
            #_, next_word = torch.max(prob, dim=1)  # for greedy sampling
            next_word = torch.multinomial(F.softmax(logits, dim=-1), 1).item()
            ys = torch.cat([ys,
                            torch.ones(1, 1).type_as(src.data).fill_(next_word)], dim=0)
            if next_word == EOS_IDX:
                break
        return ys