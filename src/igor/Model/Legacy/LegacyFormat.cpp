#include <igor/Model/Legacy/LegacyFormat.h>

namespace igor::model::legacy {

Seq_type_String legacy_gene_class_to_seq_type(Gene_class_legacy gene_class)
{
    switch (gene_class) {
    case V_gene_legacy:
        return "V_gene_seq";
    case D_gene_legacy:
        return "D_gene_seq";
    case J_gene_legacy:
        return "J_gene_seq";
    case VD_genes:
        return "VD_ins_seq";
    case DJ_genes:
        return "DJ_ins_seq";
    case VJ_genes:
        return "VJ_ins_seq";
    default:
        return "Undefined_seq";
    }
}

Gene_class_legacy legacy_gene_class_of(Event_type event_type, const Seq_type_String &seq_type)
{
    if (event_type == GeneChoice_t || event_type == Deletion_t) {
        if (seq_type == "V_gene_seq")
            return V_gene_legacy;
        if (seq_type == "D_gene_seq")
            return D_gene_legacy;
        if (seq_type == "J_gene_seq")
            return J_gene_legacy;
    } else if (event_type == Insertion_t || event_type == Dinuclmarkov_t) {
        if (seq_type == "VD_ins_seq")
            return VD_genes;
        if (seq_type == "DJ_ins_seq")
            return DJ_genes;
        if (seq_type == "VJ_ins_seq")
            return VJ_genes;
    }
    return Undefined_gene_legacy;
}

} // namespace igor::model::legacy
